import sys
from natsort import natsorted
from collections import OrderedDict
from .utils import zopen, check_inputs, which2
from .fasta import softwrap, RevComp, fasta2lengths
from .gff import gff2dict, dict2gff3, dict2gtf, gtf2dict
from .genbank import tbl2dict, dict2tbl, table2asn
import uuid
import os


def convert(args):
    check_inputs([args.input] + [args.fasta])
    if not args.input_format:  # we have to guess
        if args.input.endswith((".tbl", ".tbl.gz")):
            args.input_format = "tbl"
        elif args.input.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            args.input_format = "gff3"
        elif args.input.endswith((".gtf")):
            args.input_format = "gtf"
        else:
            sys.stderr.write(
                "Error: unable to determine -i,--input format: {}".format(args.input)
            )
            raise SystemExit(1)
    if not args.output_format:  # guess again
        if args.out.endswith((".tbl", ".tbl.gz")):
            args.output_format = "tbl"
        elif args.out.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            args.output_format = "gff3"
        elif args.out.endswith((".gtf")):
            args.output_format = "gtf"
        elif args.out.endswith((".gbk", ".gb", ".gbf", ".gbff", ".gbk")):
            args.output_format = "gbff"
            if not which2("table2asn"):
                sys.stderr.write(
                    "ERROR: table2asn is not in PATH and is required to generate genbank output"
                )
                raise SystemExit(1)
        else:
            sys.stderr.write(
                "Error: unable to determine -o,--output format: {}".format(args.out)
            )
            raise SystemExit(1)
    # okay now we can load and convert
    if args.input_format == "tbl":
        if args.output_format == "gff3":
            tbl2gff3(args.input, args.fasta, output=args.out, table=1)
        elif args.output_format == "proteins":
            tbl2proteins(
                args.input,
                args.fasta,
                output=args.out,
                table=1,
                strip_stop=args.strip_stop,
            )
        elif args.output_format == "gtf":
            tbl2gtf(args.input, args.fasta, output=args.out, table=1)
        elif args.output_format == "transcripts":
            tbl2transcripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "cds-transcripts":
            tbl2cdstranscripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "gbff":
            tbl2gbff(
                args.input,
                args.fasta,
                output=args.out,
                organism=args.organism,
                strain=args.strain,
            )
    elif args.input_format == "gff3":
        if args.output_format == "tbl":
            gff2tbl(args.input, args.fasta, output=args.out, table=1, debug=args.debug)
        elif args.output_format == "proteins":
            gff2proteins(
                args.input,
                args.fasta,
                output=args.out,
                table=1,
                strip_stop=args.strip_stop,
                debug=args.debug,
            )
        elif args.output_format == "gtf":
            gff2gtf(args.input, args.fasta, output=args.out, table=1, debug=args.debug)
        elif args.output_format == "transcripts":
            gff2transcripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "cds-transcripts":
            gff2cdstranscripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "gbff":
            gff2gbff(
                args.input,
                args.fasta,
                output=args.out,
                organism=args.organism,
                strain=args.strain,
            )
    elif args.input_format == "gtf":
        if args.output_format == "tbl":
            gtf2tbl(args.input, args.fasta, output=args.out, table=1, debug=args.debug)
        elif args.output_format == "gff3":
            gtf2gff(args.input, args.fasta, output=args.out, table=1, debug=args.debug)
        elif args.output_format == "proteins":
            gtf2proteins(
                args.input,
                args.fasta,
                output=args.out,
                table=1,
                strip_stop=args.strip_stop,
                debug=args.debug,
            )
        elif args.output_format == "transcripts":
            gtf2transcripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "cds-transcripts":
            gtf2cdstranscripts(
                args.input,
                args.fasta,
                output=args.out,
            )
        elif args.output_format == "gbff":
            gtf2gbff(
                args.input,
                args.fasta,
                output=args.out,
                organism=args.organism,
                strain=args.strain,
            )

    elif args.input_format == "miniprot":  # this is an alt GFF3 format
        pass


def _dict2proteins(input, output=False, strip_stop=False):
    if output:
        protout = zopen(output, "w")
    else:
        protout = sys.stdout
    for k, v in natsorted(list(input.items())):
        if "pseudo" in v:
            if v["pseudo"]:
                continue
        for i in range(0, len(v["ids"])):
            if v["type"][i] == "mRNA" and not v["CDS"][i]:
                continue
            if v["type"][i] == "mRNA" and not len(v["ids"]) == len(v["mRNA"]) == len(
                v["CDS"]
            ):
                continue
            if v["type"][i] == "mRNA":
                try:
                    Prot = v["protein"][i]
                except IndexError:
                    print(k, v)
                    raise SystemExit(1)
                if strip_stop:
                    Prot = Prot.rstrip("*")
                protout.write(">{:} {:}\n{:}\n".format(v["ids"][i], k, softwrap(Prot)))
    if output:
        protout.close()


def _dict2transcripts(input, output=False):
    """
    function to generate transcripts from dictionary
    """
    if output:
        tranout = zopen(output, "w")
    else:
        tranout = sys.stdout
    # write transcripts
    for k, v in natsorted(list(input.items())):
        for i, x in enumerate(v["ids"]):
            try:
                Transcript = str(v["transcript"][i])
                if v["strand"] == "-":
                    Transcript = RevComp(Transcript)
                tranout.write(">{:} {:}\n{:}\n".format(x, k, softwrap(Transcript)))
            except IndexError:
                pass
    if output:
        tranout.close()


def _dict2cdstranscripts(input, output=False):
    """
    function to generate CDS transcripts from dictionary
    """
    if output:
        tranout = zopen(output, "w")
    else:
        tranout = sys.stdout
    # write transcripts
    for k, v in natsorted(list(input.items())):
        for i, x in enumerate(v["ids"]):
            try:
                Transcript = str(v["cds_transcript"][i])
                if v["strand"] == "-":
                    Transcript = RevComp(Transcript)
                tranout.write(">{:} {:}\n{:}\n".format(x, k, softwrap(Transcript)))
            except IndexError:
                pass
    if output:
        tranout.close()


def tbl2gff3(tbl, fasta, output=False, table=1):
    """Convert NCBI TBL format to GFF3 format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write to GFF3 output. Default is to write to stdout.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        annotation file in GFF3 format

    """
    # load into dictionary
    Genes, parse_errors = tbl2dict(tbl, fasta, table=table)
    # write to GenBank format
    dict2gff3(Genes, output=output)


def tbl2gbff(tbl, fasta, output=False, table=1, organism=False, strain=False):
    """Convert NCBI TBL format to GenBank format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write to GenBank output.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        annotation file in GenBank format

    """
    # write to Genbank format via table2asn
    table2asn(tbl, fasta, output=output, table=table, organism=organism, strain=strain)


def tbl2gtf(tbl, fasta, output=False, table=1):
    """Convert NCBI TBL format to GTF format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write to GTF output.
    Only coding genes are output with this method. Default is to write to stdout.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        annotation file in GTF format

    """
    # load annotation
    Genes, parse_errors = tbl2dict(tbl, fasta, table=table)
    # write to GTF format
    dict2gtf(Genes, output=output)


def tbl2proteins(tbl, fasta, output=False, table=1, strip_stop=False):
    """Convert NCBI TBL format to translated protein FASTA format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write protein coding
    translations to FASTA format.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    strip_stop : bool, default=False
        remove stop codons (*) from translation
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load annotation
    Genes, parse_errors = tbl2dict(tbl, fasta, table=table)
    # write to protein fasta
    _dict2proteins(Genes, output=output, strip_stop=strip_stop)


def gff2proteins(gff, fasta, output=False, table=1, strip_stop=False, debug=False):
    """Convert GFF3 format to translated protein FASTA format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write protein coding
    translations to FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GFF3 format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    strip_stop : bool, default=False
        remove stop codons (*) from translation
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2proteins(Genes, output=output, strip_stop=strip_stop)


def gff2transcripts(gff, fasta, output=False, table=1, debug=False):
    """Convert GFF3 format to transcript FASTA format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write
    transcripts in FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GFF3 format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2transcripts(Genes, output=output)


def gff2cdstranscripts(gff, fasta, output=False, table=1, debug=False):
    """Convert GFF3 format to CDS transcript [no UTRs] FASTA format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write
    CDS transcripts in FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GFF3 format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2cdstranscripts(Genes, output=output)


def tbl2transcripts(tbl, fasta, output=False, table=1):
    """Convert NCBI TBL format to transcript FASTA format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write
    transcripts in FASTA format.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load annotation
    Genes, parse_errors = tbl2dict(tbl, fasta, table=table)
    # write to protein fasta
    _dict2transcripts(Genes, output=output)


def tbl2cdstranscripts(tbl, fasta, output=False, table=1):
    """Convert NCBI TBL format to CDS transcript [no UTRS] in FASTA format.

    Will parse NCBI TBL format into GFFtk annotation dictionary and then write
    CDS transcripts in FASTA format.

    Parameters
    ----------
    tbl : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load annotation
    Genes, parse_errors = tbl2dict(tbl, fasta, table=table)
    # write to protein fasta
    _dict2cdstranscripts(Genes, output=output)


def gff2gtf(gff, fasta, output=False, table=1, debug=False):
    """Convert GFF3 format to GTF format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write to GTF output.
    Only coding genes are output with this method. Default is to write to stdout.

    Parameters
    ----------
    gff : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        annotation file in GTF format

    """
    # load annotation
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # write to GTF format
    dict2gtf(Genes, output=output)


def gff2gbff(
    gff, fasta, output=False, table=1, organism=False, strain=False, debug=False
):
    """Convert GFF3 format to GenBank format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write to GenBank output.

    Parameters
    ----------
    gff : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        annotation file in GenBank format

    """
    # write to tbl format
    tmpTbl = f"{uuid.uuid4()}.tbl"

    def _sortDict(d):
        return (d[1]["location"][0], d[1]["location"][1])

    # load gff into dictionary
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # now sort the dictionary
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v["contig"] in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)
    # get contig lengths
    scaffLen = fasta2lengths(fasta)
    # finally write output
    errors, duplicates, pseudo, nocds = dict2tbl(
        sortedGenes,
        scaff2genes,
        scaffLen,
        "CFMR",
        "12345",
        [],
        output=tmpTbl,
        annotations=True,
        external=True,
    )
    # write to Genbank format via table2asn
    table2asn(
        tmpTbl, fasta, output=output, table=table, organism=organism, strain=strain
    )
    os.remove(tmpTbl)


def gtf2gff(gff, fasta, output=False, table=1, debug=False):
    """Convert GTF format to GFF format.

    Will parse GTF format into GFFtk annotation dictionary and then write to GFF3 output.
    Only coding genes are output with this method. Default is to write to stdout.

    Parameters
    ----------
    gff : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        annotation file in GTF format

    """
    # load annotation
    Genes = gtf2dict(gff, fasta, table=table, debug=debug)
    # write to GTF format
    dict2gff3(Genes, output=output)


def gtf2gbff(
    gtf, fasta, output=False, table=1, organism=False, strain=False, debug=False
):
    """Convert GTF format to GenBank format.

    Will parse GTF format into GFFtk annotation dictionary and then write to GenBank output.

    Parameters
    ----------
    gtf : filename
        genome annotation text file in NCBI tbl format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    output : str, default=sys.stdout
        annotation file in GenBank format

    """
    # write to tbl format
    tmpTbl = f"{uuid.uuid4()}.tbl"

    def _sortDict(d):
        return (d[1]["location"][0], d[1]["location"][1])

    # load gff into dictionary
    Genes = gtf2dict(gtf, fasta, table=table, debug=debug)
    # now sort the dictionary
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v["contig"] in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)
    # get contig lengths
    scaffLen = fasta2lengths(fasta)
    # finally write output
    errors, duplicates, pseudo, nocds = dict2tbl(
        sortedGenes,
        scaff2genes,
        scaffLen,
        "CFMR",
        "12345",
        [],
        output=tmpTbl,
        annotations=True,
        external=True,
    )
    # write to Genbank format via table2asn
    table2asn(
        tmpTbl, fasta, output=output, table=table, organism=organism, strain=strain
    )
    os.remove(tmpTbl)


def gff2tbl(gff, fasta, output=False, table=1, debug=False):
    """Convert GFF3 format to NCBI TBL format .

    Will parse GFF3 annotation format into GFFtk annotation dictionary and then write to
    NCBI TBL output. Default is to write to stdout.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GFF3 format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        annotation file in NCBI TBL format

    """

    def _sortDict(d):
        return (d[1]["location"][0], d[1]["location"][1])

    # load gff into dictionary
    Genes = gff2dict(gff, fasta, table=table, debug=debug)
    # now sort the dictionary
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v["contig"] in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)
    # get contig lengths
    scaffLen = fasta2lengths(fasta)
    # finally write output
    errors, duplicates, pseudo, nocds = dict2tbl(
        sortedGenes,
        scaff2genes,
        scaffLen,
        "CFMR",
        "12345",
        [],
        output=output,
        annotations=True,
        external=True,
    )


def gtf2tbl(gff, fasta, output=False, table=1, debug=False):
    """Convert GTF format to NCBI TBL format .

    Will parse GTF annotation format into GFFtk annotation dictionary and then write to
    NCBI TBL output. Default is to write to stdout.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GTF format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        annotation file in NCBI TBL format

    """

    def _sortDict(d):
        return (d[1]["location"][0], d[1]["location"][1])

    # load gff into dictionary
    Genes = gtf2dict(gff, fasta, table=table, debug=debug)
    # now sort the dictionary
    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v["contig"] in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)
    # get contig lengths
    scaffLen = fasta2lengths(fasta)
    # finally write output
    errors, duplicates, pseudo, nocds = dict2tbl(
        sortedGenes,
        scaff2genes,
        scaffLen,
        "CFMR",
        "12345",
        [],
        output=output,
        annotations=True,
        external=True,
    )


def gtf2proteins(gff, fasta, output=False, table=1, strip_stop=False, debug=False):
    """Convert GTF format to translated protein FASTA format.

    Will parse GTF format into GFFtk annotation dictionary and then write protein coding
    translations to FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GTF format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    strip_stop : bool, default=False
        remove stop codons (*) from translation
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gtf2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2proteins(Genes, output=output, strip_stop=strip_stop)


def gtf2transcripts(gff, fasta, output=False, table=1, debug=False):
    """Convert GTF format to transcript FASTA format.

    Will parse GTF format into GFFtk annotation dictionary and then write
    transcripts in FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GTF format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gtf2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2transcripts(Genes, output=output)


def gtf2cdstranscripts(gff, fasta, output=False, table=1, debug=False):
    """Convert GTF format to CDS transcript [no UTRs] FASTA format.

    Will parse GFF3 format into GFFtk annotation dictionary and then write
    CDS transcripts in FASTA format.

    Parameters
    ----------
    gff : filename
        genome annotation text file in GTF format
    fasta : filename
        genome sequence in FASTA format
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    output : str, default=sys.stdout
        translated amino acids (proteins) in FASTA format

    """
    # load gff into dictionary
    Genes = gtf2dict(gff, fasta, table=table, debug=debug)
    # write to protein fasta
    _dict2cdstranscripts(Genes, output=output)

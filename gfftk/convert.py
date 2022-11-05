import sys
from natsort import natsorted
from collections import OrderedDict
from .utils import zopen, check_inputs
from .fasta import softwrap, RevComp, fasta2lengths
from .gff import gff2dict, dict2gff3, dict2gtf
from .genbank import tbl2dict, dict2tbl


def convert(args):
    check_inputs([args.input] + [args.fasta])
    if not args.input_format:  # we have to guess
        if args.input.endswith((".tbl", ".tbl.gz")):
            args.input_format = "tbl"
        elif args.input.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            args.input_format = "gff3"
        else:
            sys.stderr.write(
                "Error: unable to determine -i,--input format: {}".format(args.input)
            )
            sys.exit(1)
    if not args.output_format:  # guess again
        if args.out.endswith((".tbl", ".tbl.gz")):
            args.output_format = "tbl"
        elif args.out.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
            args.output_format = "gff3"
        else:
            sys.stderr.write(
                "Error: unable to determine -o,--output format: {}".format(args.out)
            )
            sys.exit(1)
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
        elif args.output_format == "nucleotides":
            if not args.out:
                sys.stderr.write(
                    "Error: if output format is nucleotides, then -o,--out is required as basename for output files"
                )
                sys.exit(1)

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
        elif args.output_format == "nucleotides":
            if not args.out:
                sys.stderr.write(
                    "Error: if output format is nucleotides, then -o,--out is required as basename for output files"
                )
                sys.exit(1)
    elif args.input_format == 'miniprot':  # this is an alt GFF3 format
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
                    sys.exit(1)
                if strip_stop:
                    Prot = Prot.rstrip("*")
                protout.write(">{:} {:}\n{:}\n".format(v["ids"][i], k, softwrap(Prot)))
    if output:
        protout.close()


def _dict2nucleotides(input, prots=False, trans=False, cdstrans=False):
    """
    function to generate protein and transcripts from dictionary
    """
    if prots:
        protout = zopen(prots, mode="w")
    if trans:
        tranout = zopen(trans, mode="w")
    if cdstrans:
        cdsout = zopen(cdstrans, mode="w")
    # write to protein and transcripts
    for k, v in natsorted(list(input.items())):
        for i, x in enumerate(v["ids"]):
            if trans:
                try:
                    Transcript = str(v["transcript"][i])
                    if v["strand"] == "-":
                        Transcript = RevComp(Transcript)
                    tranout.write(">{:} {:}\n{:}\n".format(x, k, softwrap(Transcript)))
                except IndexError:
                    pass
            if cdstrans:
                try:
                    CDStranscript = str(v["cds_transcript"][i])
                    if v["strand"] == "-":
                        CDStranscript = RevComp(CDStranscript)
                    cdsout.write(
                        ">{:} {:}\n{:}\n".format(x, k, softwrap(CDStranscript))
                    )
                except IndexError:
                    pass
            if prots:
                if v["type"][i] == "mRNA":
                    try:
                        Prot = v["protein"][i]
                    except IndexError:
                        sys.stderr.write("ERROR writing protein: {}\n{}\n".format(k, v))
                        sys.exit(1)
                    Prot = Prot.rstrip("*")
                    protout.write(">{:} {:}\n{:}\n".format(x, k, softwrap(Prot)))
    if prots:
        protout.close()
    if trans:
        tranout.close()
    if cdstrans:
        cdsout.close()


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
    # write to GFF3 format
    dict2gff3(Genes, output=output)


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

#def miniprot2gff3
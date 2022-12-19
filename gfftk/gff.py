import sys
import random
import concurrent.futures
from collections import OrderedDict
from natsort import natsorted
from .fasta import translate, fasta2dict, fasta2headers, getSeqRegions
from .utils import zopen
from urllib.parse import unquote
import uuid


def start_end_gap(seq, coords):
    if seq.startswith("N"):
        oldLen = len(seq)
        seq = seq.lstrip("N")
        numLeftStripped = oldLen - len(seq)
        coords[0] = (coords[0][0] + numLeftStripped, coords[0][1])
    if seq.endswith("N"):
        oldLen = len(seq)
        seq = seq.rstrip("N")
        numRightStripped = oldLen - len(seq)
        coords[-1] = (coords[-1][0], coords[-1][1] - numRightStripped)
    return seq, coords


def _gff_default_parser(gff, fasta, Genes):
    # this is the default general parser to populate the dictionary
    # idea is to go through line by line and parse the records and add to Genes dictionary
    errors = {
        "contig_name": [],
        "columns": [],
        "comments": [],
        "unparsed_attributes": [],
        "no_parent": [],
        "no_id": [],
    }
    idParent = {}
    SeqRecords = fasta2headers(fasta)
    with zopen(gff) as input:
        for line in input:
            if line.startswith("\n") or line.startswith("#"):
                errors["comments"].append(line)
                continue
            line = line.rstrip()
            # skip lines that aren't 9 columns
            if not line.count("\t") == 8:
                errors["columns"].append(line)
                continue
            (
                contig,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = line.split("\t")
            if feature not in [
                "gene",
                "mRNA",
                "transcript",
                "exon",
                "CDS",
                "tRNA",
                "ncRNA",
                "rRNA",
                "pseudogene",
                "five_prime_UTR",
                "five_prime_utr",
                "three_prime_UTR",
                "three_prime_utr",
            ]:
                continue
            if not contig in SeqRecords:
                errors["contig_name"].append(line)
                continue
            attributes = unquote(attributes)
            source = unquote(source)
            feature = unquote(feature)
            start = int(start)
            end = int(end)
            ID, Parent, Name, Product, GeneFeature, gbkey = (None,) * 6
            info = {}
            for field in attributes.split(";"):
                try:
                    k, v = field.split("=", 1)
                    info[k] = v.strip()
                except (IndexError, ValueError) as E:
                    pass
            # now can lookup in info dict for values
            ID = info.get("ID", None)
            Parent = info.get("Parent", None)
            Name = info.get("Name", None)
            DBxref = info.get("DBxref", None)
            if not DBxref and "Dbxref" in info:
                DBxref = info.get("Dbxref")
            if not DBxref and "dbxref" in info:
                DBxref = info.get("dbxref")
            if DBxref:
                if "," in DBxref:
                    DBxref = DBxref.split(",")
                else:
                    DBxref = [DBxref]
            else:
                DBxref = []
            GO = info.get("Ontology_term", None)
            if GO:
                if "," in GO:
                    GO = GO.split(",")
                else:
                    GO = [GO]
            else:
                GO = []
            ECnum = info.get("EC_number", None)
            if ECnum:
                if "," in ECnum:
                    ECnum = ECnum.split(",")
                else:
                    ECnum = [ECnum]
            else:
                ECnum = []
            Note = info.get("Note", None)
            if not Note and "note" in info:
                Note = info.get("note")
            if Note:
                if "," in Note:
                    Note = Note.split(",")
                else:
                    Note = [Note]
            else:
                Note = []
            Product = info.get("Product", None)
            if not Product and "product" in info:
                Product = info.get("product")
            if not Product and "description" in info:
                Product = info.get("description")
            synonyms = info.get("Alias", None)
            if not synonyms and "gene_synonym" in info:
                synonyms = info.get("gene_synonym")
            if synonyms:
                if "," in synonyms:
                    synonyms = synonyms.split(",")
                else:
                    synonyms = [synonyms]
            else:
                synonyms = []
            gbkey = info.get("gbkey", None)
            # for error reporting capture unparsed keys
            for attr, value in info.items():
                if not attr in [
                    "ID",
                    "Parent",
                    "Name",
                    "DBxref",
                    "Dbxref",
                    "dbxref",
                    "Ontology_term",
                    "EC_number",
                    "Note",
                    "note",
                    "Product",
                    "product",
                    "description",
                    "Alias",
                    "gbkey",
                    "gene_synonym",
                ]:
                    if not attr in errors["unparsed_attributes"]:
                        errors["unparsed_attributes"].append(attr)
            # now we can do add to dictionary these parsed values
            # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
            if feature in ["gene", "pseudogene"]:
                if not ID in Genes:
                    if feature == "pseudogene":
                        pseudoFlag = True
                    else:
                        pseudoFlag = False
                    Genes[ID] = {
                        "name": Name,
                        "type": [],
                        "transcript": [],
                        "cds_transcript": [],
                        "protein": [],
                        "5UTR": [],
                        "3UTR": [],
                        "gene_synonym": synonyms,
                        "codon_start": [],
                        "ids": [],
                        "CDS": [],
                        "mRNA": [],
                        "strand": strand,
                        "EC_number": [],
                        "location": (start, end),
                        "contig": contig,
                        "product": [],
                        "source": source,
                        "phase": [],
                        "db_xref": [],
                        "go_terms": [],
                        "note": [],
                        "partialStart": [],
                        "partialStop": [],
                        "pseudo": pseudoFlag,
                    }
                else:
                    if start < Genes[ID]["location"][0]:
                        Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                    if end > Genes[ID]["location"][1]:
                        Genes[ID]["location"] = (Genes[ID]["location"][0], end)
            else:
                if not ID:
                    if Name:
                        ID = Name
                    # one of the dumbest things I've seen in ensembl they have only Parent=
                    elif feature in ["three_prime_UTR", "five_prime_UTR"]:
                        ID = str(uuid.uuid4())
                    else:
                        errors["no_id"].append(line)
                        continue
                if not Parent:
                    errors["no_parent"].append(line)
                    continue
                if feature in ["mRNA", "transcript", "tRNA", "ncRNA", "rRNA"]:
                    if gbkey and gbkey == "misc_RNA":
                        feature = "ncRNA"
                    if not Product:
                        if feature in ["mRNA", "transcript"]:
                            Product = "hypothetical protein"
                    if not Parent in Genes:
                        Genes[Parent] = {
                            "name": Name,
                            "type": [feature],
                            "transcript": [],
                            "cds_transcript": [],
                            "protein": [],
                            "5UTR": [[]],
                            "3UTR": [[]],
                            "codon_start": [[]],
                            "ids": [ID],
                            "CDS": [[]],
                            "mRNA": [[]],
                            "strand": strand,
                            "location": (start, end),
                            "contig": contig,
                            "product": [Product],
                            "source": source,
                            "phase": [[]],
                            "gene_synonym": synonyms,
                            "db_xref": [DBxref],
                            "go_terms": [GO],
                            "EC_number": [ECnum],
                            "note": [Note],
                            "partialStart": [False],
                            "partialStop": [False],
                            "pseudo": False,
                        }
                    else:
                        Genes[Parent]["ids"].append(ID)
                        Genes[Parent]["mRNA"].append([])
                        Genes[Parent]["CDS"].append([])
                        Genes[Parent]["phase"].append([])
                        Genes[Parent]["5UTR"].append([])
                        Genes[Parent]["3UTR"].append([])
                        Genes[Parent]["codon_start"].append([])
                        Genes[Parent]["partialStart"].append(False)
                        Genes[Parent]["partialStop"].append(False)
                        Genes[Parent]["product"].append(Product)
                        Genes[Parent]["db_xref"].append(DBxref)
                        Genes[Parent]["EC_number"].append(ECnum)
                        Genes[Parent]["gene_synonym"] += synonyms
                        Genes[Parent]["go_terms"].append(GO)
                        Genes[Parent]["note"].append(Note)
                        Genes[Parent]["type"].append(feature)
                        # double check mRNA features are contained in gene coordinates
                        if start < Genes[Parent]["location"][0]:
                            Genes[Parent]["location"] = (
                                start,
                                Genes[Parent]["location"][1],
                            )
                        if end > Genes[Parent]["location"][1]:
                            Genes[Parent]["location"] = (
                                Genes[Parent]["location"][0],
                                end,
                            )
                    if not ID in idParent:
                        idParent[ID] = Parent
                # treat exon features
                elif feature == "exon":
                    if "," in Parent:
                        parents = Parent.split(",")
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {
                                    "name": Name,
                                    "type": [],
                                    "transcript": [],
                                    "cds_transcript": [],
                                    "protein": [],
                                    "5UTR": [[]],
                                    "3UTR": [[]],
                                    "codon_start": [[]],
                                    "ids": [p],
                                    "CDS": [],
                                    "mRNA": [[(start, end)]],
                                    "strand": strand,
                                    "location": None,
                                    "contig": contig,
                                    "product": [],
                                    "source": source,
                                    "phase": [[]],
                                    "db_xref": [],
                                    "go_terms": [],
                                    "EC_number": [],
                                    "note": [],
                                    "partialStart": [False],
                                    "partialStop": [False],
                                    "pseudo": False,
                                    "gene_synonym": synonyms,
                                }
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]["ids"].index(p)
                                Genes[GeneFeature]["mRNA"][i].append((start, end))
                # treat codings sequence features
                elif feature == "CDS":
                    if "," in Parent:
                        parents = Parent.split(",")
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {
                                    "name": Name,
                                    "type": [],
                                    "transcript": [],
                                    "cds_transcript": [],
                                    "protein": [],
                                    "5UTR": [[]],
                                    "3UTR": [[]],
                                    "codon_start": [[]],
                                    "ids": [p],
                                    "CDS": [[(start, end)]],
                                    "mRNA": [],
                                    "strand": strand,
                                    "location": None,
                                    "contig": contig,
                                    "product": [],
                                    "source": source,
                                    "phase": [[]],
                                    "db_xref": [],
                                    "go_terms": [],
                                    "EC_number": [],
                                    "note": [],
                                    "partialStart": [False],
                                    "partialStop": [False],
                                    "pseudo": False,
                                    "gene_synonym": synonyms,
                                }
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]["ids"].index(p)
                                Genes[GeneFeature]["CDS"][i].append((start, end))
                                # add phase
                                try:
                                    Genes[GeneFeature]["phase"][i].append(int(phase))
                                except ValueError:
                                    Genes[GeneFeature]["phase"][i].append("?")
                # treat 5' UTRs
                elif feature == "five_prime_UTR" or feature == "five_prime_utr":
                    if "," in Parent:
                        parents = Parent.split(",")
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {
                                    "name": Name,
                                    "type": [],
                                    "transcript": [],
                                    "cds_transcript": [],
                                    "protein": [],
                                    "5UTR": [[(start, end)]],
                                    "3UTR": [[]],
                                    "codon_start": [[]],
                                    "ids": [p],
                                    "CDS": [],
                                    "mRNA": [[(start, end)]],
                                    "strand": strand,
                                    "location": None,
                                    "contig": contig,
                                    "product": [],
                                    "source": source,
                                    "phase": [[]],
                                    "db_xref": [],
                                    "go_terms": [],
                                    "EC_number": [],
                                    "note": [],
                                    "partialStart": [False],
                                    "partialStop": [False],
                                    "pseudo": False,
                                    "gene_synonym": synonyms,
                                }
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]["ids"].index(p)
                                Genes[GeneFeature]["5UTR"][i].append((start, end))
                # treat 3' UTR
                elif feature == "three_prime_UTR" or feature == "three_prime_utr":
                    if "," in Parent:
                        parents = Parent.split(",")
                    else:
                        parents = [Parent]
                    for p in parents:
                        if p in idParent:
                            GeneFeature = idParent.get(p)
                        if GeneFeature:
                            if not GeneFeature in Genes:
                                Genes[GeneFeature] = {
                                    "name": Name,
                                    "type": [],
                                    "transcript": [],
                                    "cds_transcript": [],
                                    "protein": [],
                                    "5UTR": [[]],
                                    "3UTR": [[(start, end)]],
                                    "codon_start": [[]],
                                    "ids": [p],
                                    "CDS": [],
                                    "mRNA": [[(start, end)]],
                                    "strand": strand,
                                    "location": None,
                                    "contig": contig,
                                    "product": [],
                                    "source": source,
                                    "phase": [[]],
                                    "db_xref": [],
                                    "go_terms": [],
                                    "EC_number": [],
                                    "note": [],
                                    "partialStart": [False],
                                    "partialStop": [False],
                                    "pseudo": False,
                                    "gene_synonym": synonyms,
                                }
                            else:
                                # determine which transcript this is get index from id
                                i = Genes[GeneFeature]["ids"].index(p)
                                Genes[GeneFeature]["3UTR"][i].append((start, end))
    return Genes, errors


def validate_and_translate_models(
    k, v, SeqRecords, gap_filter=False, table=1, logger=sys.stderr.write
):
    # take the Genes dictionary and validate gene models
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    # return sorted mRNA, sorted CDS, transcripts, translations, proper phase
    results = {
        "gene_synonym": [],
        "location": v["location"],
        "mRNA": [],
        "CDS": [],
        "protein": [],
        "transcript": [],
        "cds_transcript": [],
        "codon_start": [],
        "partialStart": [],
        "partialStop": [],
        "type": [],
    }
    assert len(v["ids"]) == len(v["type"])
    for i in range(0, len(v["ids"])):
        if v["type"][i] in ["mRNA", "tRNA", "ncRNA", "rRNA", "transcript"]:
            if v["strand"] == "+":
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0], reverse=True)
            mrnaSeq = getSeqRegions(SeqRecords, v["contig"], sortedExons)
            if gap_filter:
                mrnaSeq, sortedExons = start_end_gap(mrnaSeq, sortedExons)
            results["mRNA"].append(sortedExons)
            results["transcript"].append(mrnaSeq)
        if v["type"][i] in ["mRNA", "transcript"]:
            if not v["CDS"][i]:
                results["CDS"].append([])
                results["type"].append("ncRNA")
                results["protein"].append(None)
                results["cds_transcript"].append(None)
                results["codon_start"].append(None)
                results["partialStart"].append(None)
                results["partialStop"].append(None)
            else:
                results["type"].append("mRNA")
                if v["strand"] == "+":
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(
                        v["CDS"][i], key=lambda tup: tup[0], reverse=True
                    )
                # get the codon_start by getting first CDS phase + 1
                indexStart = [
                    x for x, y in enumerate(v["CDS"][i]) if y[0] == sortedCDS[0][0]
                ]
                cdsSeq = getSeqRegions(SeqRecords, v["contig"], sortedCDS)
                if gap_filter:
                    cdsSeq, v["CDS"][i] = start_end_gap(cdsSeq, v["CDS"][i])
                protSeq, codon_start = (None,) * 2
                if (
                    "?" in v["phase"][i]
                ):  # dont know the phase -- malformed GFF3, try to find best CDS
                    translateResults = []
                    for y in [1, 2, 3]:
                        protSeq = translate(cdsSeq, v["strand"], y - 1, table=table)
                        if not protSeq:
                            logger(
                                "Translation of {} using {} phase failed\n".format(
                                    v["ids"][i], y - 1
                                )
                            )
                            sys.exit(1)
                        numStops = protSeq.count("*")
                        if protSeq[-1] == "*":
                            numStops -= 1
                        translateResults.append((y, numStops, protSeq))
                    sortedResults = sorted(translateResults, key=lambda tup: tup[1])
                    codon_start = sortedResults[0][0]
                    protSeq = sortedResults[0][2]
                else:
                    try:
                        codon_start = int(v["phase"][i][indexStart[0]]) + 1
                    except IndexError:
                        pass
                    # translate and get protein sequence
                    protSeq = translate(
                        cdsSeq, v["strand"], codon_start - 1, table=table
                    )
                results["codon_start"].append(codon_start)
                if codon_start > 1:
                    if v["strand"] == "+":
                        cdsSeq = cdsSeq[codon_start - 1 :]
                    elif v["strand"] == "-":
                        endTrunc = len(cdsSeq) - codon_start - 1
                        cdsSeq = cdsSeq[0:endTrunc]
                    else:
                        logger(
                            "ERROR nonsensical strand ({}) for gene {}\n".format(
                                [v["strand"], v["ids"][i]]
                            )
                        )
                results["cds_transcript"].append(cdsSeq)
                results["CDS"].append(sortedCDS)
                results["protein"].append(protSeq)
                if protSeq:
                    if protSeq.endswith("*"):
                        results["partialStop"].append(False)
                    else:
                        results["partialStop"].append(True)
                    if codon_start == 1 and protSeq.startswith("M"):
                        results["partialStart"].append(False)
                    else:
                        results["partialStart"].append(True)
        else:
            results["CDS"].append([])
            results["type"].append(v["type"][i])
            results["codon_start"].append(None)
            results["partialStart"].append(None)
            results["partialStop"].append(None)
            results["protein"].append(None)
            results["cds_transcript"].append(None)
        # since its possible updated the mRNA/CDS fields, double check that gene coordinates are ok
        all_mRNA_coords = [item for sublist in results["mRNA"] for item in sublist]
        try:
            results["location"] = (
                min(all_mRNA_coords, key=lambda item: item[0])[0],
                max(all_mRNA_coords, key=lambda item: item[1])[1],
            )
        except ValueError:
            continue
        # clean up any repeated synonym
        if len(v["gene_synonym"]) > 1:
            uniqueSynonyms = set(v["gene_synonym"])
            results["gene_synonym"] = list(uniqueSynonyms)
    return (k, results)


def gff2dict(
    gff,
    fasta,
    annotation=False,
    table=1,
    debug=False,
    gap_filter=False,
    logger=sys.stderr.write,
):
    """Convert GFF3 and FASTA to standardized GFFtk dictionary format.

    Annotation file in GFF3 format and genome FASTA file are parsed. The result is a dictionary that is keyed
    by locus_tag (gene name) and the value is a nested dictionary containing feature information.

    Parameters
    ----------
    gff : filename : str
        annotation text file in GFF3 format
    fasta : filename : str
        genome text file in FASTA format
    annotation : dict of str
        existing annotation dictionary
    table : int, default=1
        codon table [1]
    debug : bool, default=False
        print debug information to stderr
    gap_filter : bool, default=False
        remove gene models that span gaps in sequence
    logger : handle, default=sys.stderr.write
        where to log messages to

    Returns
    -------
    annotation : dict of dict
        standardized annotation dictionary (OrderedDict) keyed by locus_tag

    """
    # some logic here to predict format eventually
    if not annotation:
        annotation = {}
    annotation, parse_errors = _gff_default_parser(gff, fasta, annotation)
    if debug:
        for err, err_v in parse_errors.items():
            if len(err_v) > 0:  # then tell user
                if err == "comments":
                    continue
                if err == "unparsed_attributes":
                    print_errors = list(err_v)
                    logger(
                        "Found {} attribute keys that were not parsed: {}\n".format(
                            len(err_v), print_errors
                        )
                    )
                else:
                    logger(
                        "Found {} errors in {}: showing 10 randomly generated lines/records\n".format(
                            len(err_v), err
                        )
                    )
                    print_errors = random.sample(err_v, 10)
                    logger("{}\n".format("\n".join(print_errors)))
    SeqRecords = fasta2dict(fasta)
    if debug:
        logger(
            "Parsed {} contigs containing {} gene models\n".format(
                len(SeqRecords), len(annotation)
            )
        )
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    # try to use multithreading here, not sure its necessary but new syntax for me with concurrent
    results = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for k, v in list(annotation.items()):
            results.append(
                executor.submit(
                    validate_and_translate_models,
                    k,
                    v,
                    SeqRecords,
                    {"gap_filter": gap_filter, "table": table, "logger": logger},
                )
            )
            # check, update = validate_and_translate_models(v, SeqRecords, gap_filter=gap_filter, table=table, logger=logger)
    # pass the updates back to the main dictionary
    for r in results:
        gene, update = r.result()
        for key, value in update.items():
            annotation[gene][key] = value
        # some assertion statements here to ensure parsing is correct
        assert_lengths_fail = []
        for z in [
            "type",
            "mRNA",
            "CDS",
            "codon_start",
            "phase",
            "5UTR",
            "3UTR",
            "protein",
            "transcript",
            "cds_transcript",
            "partialStart",
            "partialStop",
            "product",
        ]:
            if len(annotation[gene]["ids"]) != len(annotation[gene][z]):
                assert_lengths_fail.append((z, annotation[gene][z], len(annotation[gene][z])))
        if len(assert_lengths_fail) > 0:
            logger(
                "ERROR in parsing gene {}\n{}\n{}\n".format(
                    gene, assert_lengths_fail, annotation[gene]
                )
            )
            sys.exit(1)
    return annotation


def simplifyGO(inputList):
    simple = []
    for x in inputList:
        if x.startswith("GO:"):
            simple.append(x.strip())
        elif " " in x:
            simple.append(x.split(" ")[1])
    return simple


def dict2gff3(input, output=False, debug=False):
    """Convert GFFtk standardized annotation dictionary to GFF3 file.

    Annotation dictionary generated by gff2dict or tbl2dict passed as input. This function then write to GFF3 format

    Parameters
    ----------
    input : dict of dict
        standardized annotation dictionary keyed by locus_tag
    output : str, default=sys.stdout
        annotation file in GFF3 format
    debug : bool, default=False
        print debug information to stderr

    """

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    if output:
        gffout = zopen(output, mode="w")
    else:
        gffout = sys.stdout
    gffout.write("##gff-version 3\n")
    for k, v in list(sortedGenes.items()):
        if v["name"]:
            if "gene_synonym" in v and len(v["gene_synonym"]) > 0:
                gffout.write(
                    "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};Alias={:};\n".format(
                        v["contig"],
                        v["source"],
                        v["location"][0],
                        v["location"][1],
                        v["strand"],
                        k,
                        v["name"],
                        ",".join(v["gene_synonym"]),
                    )
                )
            else:
                gffout.write(
                    "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                        v["contig"],
                        v["source"],
                        v["location"][0],
                        v["location"][1],
                        v["strand"],
                        k,
                        v["name"],
                    )
                )
        else:
            if "gene_synonym" in v and len(v["gene_synonym"]) > 0:
                gffout.write(
                    "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};Alias={:};\n".format(
                        v["contig"],
                        v["source"],
                        v["location"][0],
                        v["location"][1],
                        v["strand"],
                        k,
                        ",".join(v["gene_synonym"]),
                    )
                )
            else:
                gffout.write(
                    "{:}\t{:}\tgene\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                        v["contig"],
                        v["source"],
                        v["location"][0],
                        v["location"][1],
                        v["strand"],
                        k,
                    )
                )

        for i in range(0, len(v["ids"])):
            # make sure coordinates are sorted
            if v["strand"] == "+":
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0])
                if v['type'][i] == 'mRNA':
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0])
                if "5UTR" in v and v["5UTR"][i]:
                    sortedFive = sorted(v["5UTR"][i], key=lambda tup: tup[0])
                if "3UTR" in v and v["3UTR"][i]:
                    sortedThree = sorted(v["3UTR"][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0], reverse=True)
                if v['type'][i] == 'mRNA':
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0], reverse=True)
                if "5UTR" in v and v["5UTR"][i]:
                    sortedFive = sorted(
                        v["5UTR"][i], key=lambda tup: tup[0], reverse=True
                    )
                if "3UTR" in v and v["3UTR"][i]:
                    sortedThree = sorted(
                        v["3UTR"][i], key=lambda tup: tup[0], reverse=True
                    )
            # build extra annotations for each transcript if applicable
            extraAnnotations = ""
            if "gene_synonym" in v and len(v["gene_synonym"]) > 0:
                extraAnnotations = extraAnnotations + "Alias={:};".format(
                    ",".join(v["gene_synonym"])
                )
            if len(v["go_terms"][i]) > 0:
                go_annotations = simplifyGO(v["go_terms"][i])
                extraAnnotations = extraAnnotations + "Ontology_term={:};".format(
                    ",".join(go_annotations)
                )
            if len(v["db_xref"][i]) > 0:
                extraAnnotations = extraAnnotations + "Dbxref={:};".format(
                    ",".join(v["db_xref"][i])
                )
            if "EC_number" in v and len(v["EC_number"][i]) > 0:
                extraAnnotations = extraAnnotations + "EC_number={:};".format(
                    ",".join(v["EC_number"][i])
                )
            if len(v["note"][i]) > 0:
                CleanedNote = (
                    []
                )  # need to make sure no commas or semi-colons in these data else will cause problems in parsing GFF3 output downstream
                for x in v["note"][i]:
                    if ";" in x:
                        x = x.replace(";", ".")
                    if ":" in x:
                        base, values = x.split(":", 1)
                        if not "," in values:
                            CleanedNote.append(base + ":" + values)
                        else:
                            for y in values.split(","):
                                CleanedNote.append(base + ":" + y)
                    else:
                        CleanedNote.append(x.replace(",", ""))
                extraAnnotations = extraAnnotations + "Note={:};".format(
                    ",".join(CleanedNote)
                )
            # now write mRNA feature
            gffout.write(
                "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};{:}\n".format(
                    v["contig"],
                    v["source"],
                    v["type"][i],
                    v["location"][0],
                    v["location"][1],
                    v["strand"],
                    v["ids"][i],
                    k,
                    v["product"][i],
                    extraAnnotations,
                )
            )
            if v["type"][i] in ["mRNA", "tRNA", "ncRNA"]:
                if "5UTR" in v and v["5UTR"][i]:
                    # if 5'UTR then write those first
                    num_5utrs = len(v["5UTR"][i])
                    if num_5utrs > 0:
                        for z in range(0, num_5utrs):
                            u_num = z + 1
                            gffout.write(
                                "{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                    v["contig"],
                                    v["source"],
                                    sortedFive[z][0],
                                    sortedFive[z][1],
                                    v["strand"],
                                    v["ids"][i],
                                    u_num,
                                    v["ids"][i],
                                )
                            )
                # write the exons
                num_exons = len(v["mRNA"][i])
                for x in range(0, num_exons):
                    ex_num = x + 1
                    gffout.write(
                        "{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:};\n".format(
                            v["contig"],
                            v["source"],
                            sortedExons[x][0],
                            sortedExons[x][1],
                            v["strand"],
                            v["ids"][i],
                            ex_num,
                            v["ids"][i],
                        )
                    )
                # if 3'UTR then write
                if "3UTR" in v and v["3UTR"][i]:
                    num_3utrs = len(v["3UTR"][i])
                    if num_3utrs > 0:
                        for z in range(0, num_3utrs):
                            u_num = z + 1
                            gffout.write(
                                "{:}\t{:}\tthree_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr3p{:};Parent={:};\n".format(
                                    v["contig"],
                                    v["source"],
                                    sortedThree[z][0],
                                    sortedThree[z][1],
                                    v["strand"],
                                    v["ids"][i],
                                    u_num,
                                    v["ids"][i],
                                )
                            )
            if v["type"][i] == "mRNA":
                num_cds = len(v["CDS"][i])
                # GFF3 phase is 1 less than flat file
                current_phase = v["codon_start"][i] - 1
                for y in range(0, num_cds):
                    gffout.write(
                        "{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:};\n".format(
                            v["contig"],
                            v["source"],
                            sortedCDS[y][0],
                            sortedCDS[y][1],
                            v["strand"],
                            current_phase,
                            v["ids"][i],
                            v["ids"][i],
                        )
                    )
                    current_phase = (
                        current_phase
                        - (int(sortedCDS[y][1]) - int(sortedCDS[y][0]) + 1)
                    ) % 3
                    if current_phase == 3:
                        current_phase = 0
    if output:
        gffout.close()


def dict2gtf(input, output=False):
    """Convert GFFtk standardized annotation dictionary to GTF file.

    Annotation dictionary generated by gff2dict or tbl2dict passed as input. This function
    then write to GTF format, notably this function only writes protein coding CDS features.

    Parameters
    ----------
    input : dict of dict
        standardized annotation dictionary keyed by locus_tag
    output : str, default=sys.stdout
        annotation in GTF format

    """

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    if output:
        gtfout = zopen(output, mode="w")
    else:
        gtfout = sys.stdout
    for k, v in list(sortedGenes.items()):
        for i in range(0, len(v["ids"])):
            if len(v["CDS"][i]) < 1:
                continue
            # create attributes string
            attributes = 'gene_id "{:}"; transcript_id "{:}";'.format(k, v["ids"][i])
            if len(v["5UTR"][i]) > 0:
                for utr in v["5UTR"][i]:
                    gtfout.write(
                        "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                            v["contig"],
                            v["source"],
                            "5UTR",
                            utr[0],
                            utr[1],
                            0,
                            v["strand"],
                            0,
                            attributes,
                        )
                    )
            if not v["partialStart"][i]:
                if v["strand"] == "+":
                    startCodon = (v["CDS"][i][0][0], v["CDS"][i][0][0] + 2)
                else:
                    startCodon = (v["CDS"][i][0][1] - 2, v["CDS"][i][0][1])
                gtfout.write(
                    "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                        v["contig"],
                        v["source"],
                        "start_codon",
                        startCodon[0],
                        startCodon[1],
                        0,
                        v["strand"],
                        0,
                        attributes,
                    )
                )
            current_phase = v["codon_start"][i] - 1
            for x, cds in enumerate(v["CDS"][i]):
                if v["partialStop"][
                    i
                ]:  # then just write the whole CDS as no reason to move codon back
                    gtfout.write(
                        "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                            v["contig"],
                            v["source"],
                            "CDS",
                            cds[0],
                            cds[1],
                            0,
                            v["strand"],
                            current_phase,
                            attributes,
                        )
                    )
                else:
                    current_phase = (
                        current_phase - (int(cds[1]) - int(cds[0]) + 1)
                    ) % 3
                    if current_phase == 3:
                        current_phase = 0
                    if v["strand"] == "+":
                        if x == len(v["CDS"][i]) - 1:  # this is last one
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    v["source"],
                                    "CDS",
                                    cds[0],
                                    cds[1] - 3,
                                    0,
                                    v["strand"],
                                    current_phase,
                                    attributes,
                                )
                            )
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    v["source"],
                                    "stop_codon",
                                    cds[1] - 2,
                                    cds[1],
                                    0,
                                    v["strand"],
                                    0,
                                    attributes,
                                )
                            )
                        else:
                            try:
                                gtfout.write(
                                    "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                        v["contig"],
                                        v["source"],
                                        "CDS",
                                        cds[0],
                                        cds[1],
                                        0,
                                        v["strand"],
                                        current_phase,
                                        attributes,
                                    )
                                )
                            except IndexError:
                                print(k, v)
                                sys.exit(1)
                    else:
                        if x == len(v["CDS"][i]) - 1:  # this is last one
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    v["source"],
                                    "CDS",
                                    cds[0] + 3,
                                    cds[1],
                                    0,
                                    v["strand"],
                                    current_phase,
                                    attributes,
                                )
                            )
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    v["source"],
                                    "stop_codon",
                                    cds[0],
                                    cds[0] + 2,
                                    0,
                                    v["strand"],
                                    0,
                                    attributes,
                                )
                            )
                        else:
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    v["source"],
                                    "CDS",
                                    cds[0],
                                    cds[1],
                                    0,
                                    v["strand"],
                                    current_phase,
                                    attributes,
                                )
                            )
            if len(v["3UTR"][i]) > 0:
                for utr in v["3UTR"][i]:
                    gtfout.write(
                        "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                            v["contig"],
                            v["source"],
                            "3UTR",
                            utr[0],
                            utr[1],
                            0,
                            v["strand"],
                            0,
                            attributes,
                        )
                    )
            if i == len(v["ids"]) - 1:
                gtfout.write("\n")
    if output:
        gtfout.close()

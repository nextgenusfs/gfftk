import sys
import numpy as np
from natsort import natsorted
from .fasta import getSeqRegions, fasta2dict, translate
from .utils import zopen, readBlocks2
from .go import go_term_dict
from .interlap import InterLap


def tbl2dict(input, fasta, annotation=False, table=1, debug=False):
    """
    need a method to convert directly from NCBI tbl format to several output formats
    to avoid conversion problems with GBK files that have mutliple transcripts
    if can load funannotate dictionary directly from tbl format, then can write the other
    formats directly
    """
    if not annotation:
        annotation = {}
    errors = []
    with zopen(input) as infile:
        contig = ""
        for item in readBlocks2(infile, ">Feature", "\tgene\n"):
            if item[0].startswith(">Feature"):  # this will be contig header block
                contig = item[0].rstrip().split(" ")[-1]
            else:  # these are all gene model blocks
                (
                    geneID,
                    Name,
                    type,
                    start,
                    end,
                    fivepartial,
                    threepartial,
                    strand,
                    location,
                ) = (None,) * 9
                codon_start = []
                transcriptID = []
                proteinID = []
                synonyms = []
                product = []
                phase = "?"
                first, firstpartial, second, secondpartial = (False,) * 4
                position = None
                # check number of transcripts
                tNum = 0
                for z in item:
                    if z.startswith("\t\t\ttranscript_id"):
                        tNum += 1
                if tNum > 0:
                    tNum = int(tNum / 2)
                if tNum == 0:
                    tNum = 1
                # setup lists for transcripts
                mRNA = [[] for y in range(tNum)]
                CDS = [[] for y in range(tNum)]
                note = [[] for y in range(tNum)]
                dbxref = [[] for y in range(tNum)]
                ECnum = [[] for y in range(tNum)]
                go_terms = [[] for y in range(tNum)]
                fivepartial = [
                    False,
                ] * tNum
                threepartial = [
                    False,
                ] * tNum
                currentNum = 0
                for x in item:
                    exonF, exonR, cdsF, cdsR, cols = (None,) * 5
                    if x.endswith("\tgene\n") and not position:
                        cols = x.strip().split("\t")
                        position = "gene"
                        if cols[0].startswith("<"):
                            first = int(cols[0].split("<")[-1])
                        else:
                            first = int(cols[0])
                        if cols[1].startswith(">"):
                            second = int(cols[1].split(">")[-1])
                        else:
                            second = int(cols[1])
                        if first < second:
                            start = first
                            end = second
                            strand = "+"
                        else:
                            start = second
                            end = first
                            strand = "-"
                        location = (start, end)
                    elif x.startswith("\t\t\tgene\t"):
                        Name = x.strip().split("\t")[-1]
                    elif x.startswith("\t\t\tlocus_tag\t"):
                        geneID = x.strip().split("\t")[-1]
                    elif (
                        x.endswith("\ttRNA\n")
                        and x.count("\t") == 2
                        and position == "gene"
                    ):
                        type = "tRNA"
                        position = "tRNA"
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif (
                        x.endswith("\tncRNA\n")
                        and x.count("\t") == 2
                        and position == "gene"
                    ):
                        type = "ncRNA"
                        position = "ncRNA"
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif (
                        x.endswith("\trRNA\n")
                        and x.count("\t") == 2
                        and position == "gene"
                    ):
                        type = "rRNA"
                        position = "rRNA"
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith("\tmRNA\n") and x.count("\t") == 2:
                        if position == "CDS":
                            currentNum += 1
                        elif position == "gene":
                            type = "mRNA"
                        position = "mRNA"
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif x.endswith("\tCDS\n") and x.count("\t") == 2:
                        position = "CDS"
                        cols = x.strip().split("\t")
                        cdsF = int(cols[0].replace("<", ""))
                        cdsR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            CDS[currentNum].append((cdsF, cdsR))
                        else:
                            CDS[currentNum].append((cdsR, cdsF))
                    elif x.startswith("\t\t\tcodon_start\t"):
                        cNum = int(x.strip().split("\t")[-1])
                        codon_start.append(cNum)
                        phase = cNum - 1
                    elif x.startswith("\t\t\tproduct\t") and position != "mRNA":
                        product.append(x.strip().split("\t")[-1])
                    elif x.startswith("\t\t\ttranscript_id\t"):
                        tID = x.strip().split("|")[-1]
                        if "_mrna" in tID:
                            tID = tID.replace("_mrna", "")
                        if not tID in transcriptID:
                            transcriptID.append(tID)
                    elif x.startswith("\t\t\tprotein_id\t"):
                        pID = x.strip().split("|")[-1]
                        if not pID in proteinID:
                            proteinID.append(pID)
                    elif x.startswith("\t\t\tgene_synonym\t"):
                        synonyms.append(x.strip().split("\t")[-1])
                    elif x.startswith("\t\t\tgo_"):  # go terms
                        go_terms[currentNum].append(
                            "GO:{:}".format(x.strip().split("|")[1])
                        )
                    elif x.startswith("\t\t\tnote\t"):
                        note[currentNum].append(x.strip().split("\t")[-1])
                    elif x.startswith("\t\t\tdb_xref\t"):
                        dbxref[currentNum].append(x.strip().split("\t")[-1])
                    elif x.startswith("\t\t\tEC_number\t"):
                        ECnum[currentNum].append(x.strip().split("\t")[-1])
                    elif position == "mRNA" and x.count("\t") == 1:
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif position in ["tRNA", "ncRNA", "rRNA"] and x.count("\t") == 1:
                        cols = x.strip().split("\t")
                        exonF = int(cols[0].replace("<", ""))
                        exonR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            mRNA[currentNum].append((exonF, exonR))
                        else:
                            mRNA[currentNum].append((exonR, exonF))
                    elif position == "CDS" and x.count("\t") == 1:
                        cols = x.strip().split("\t")
                        cdsF = int(cols[0].replace("<", ""))
                        cdsR = int(cols[1].replace(">", ""))
                        if strand == "+":
                            CDS[currentNum].append((cdsF, cdsR))
                        else:
                            CDS[currentNum].append((cdsR, cdsF))
                if not geneID in annotation:
                    if type in ["tRNA", "ncRNA", "rRNA"]:
                        annotation[geneID] = {
                            "name": Name,
                            "type": [type],
                            "transcript": [],
                            "cds_transcript": [],
                            "protein": [],
                            "5UTR": [[]],
                            "3UTR": [[]],
                            "codon_start": codon_start,
                            "ids": [geneID + "-T1"],
                            "CDS": CDS,
                            "mRNA": mRNA,
                            "strand": strand,
                            "gene_synonym": synonyms,
                            "location": location,
                            "contig": contig,
                            "product": product,
                            "source": "GFFtk",
                            "phase": [phase],
                            "db_xref": dbxref,
                            "go_terms": go_terms,
                            "EC_number": ECnum,
                            "note": note,
                            "partialStart": [True],
                            "partialStop": [True],
                            "pseudo": False,
                        }
                    else:
                        annotation[geneID] = {
                            "name": Name,
                            "type": [type],
                            "transcript": [],
                            "cds_transcript": [],
                            "protein": [],
                            "5UTR": [],
                            "3UTR": [],
                            "codon_start": codon_start,
                            "ids": proteinID,
                            "CDS": CDS,
                            "mRNA": mRNA,
                            "strand": strand,
                            "gene_synonym": synonyms,
                            "location": location,
                            "contig": contig,
                            "product": product,
                            "source": "GFFtk",
                            "phase": [phase],
                            "db_xref": dbxref,
                            "go_terms": go_terms,
                            "EC_number": ECnum,
                            "note": note,
                            "partialStart": fivepartial,
                            "partialStop": threepartial,
                            "pseudo": False,
                        }
    # now we need to sort coordinates, get protein/transcript sequences and capture UTRs
    SeqRecords = fasta2dict(fasta)
    for k, v in list(annotation.items()):
        # @nextgenusfs we should clarify or rename this variable to indicate
        # i is the i-th transcript, right??
        for i in range(0, len(v["ids"])):
            print(v['type'], i)
            if v["type"][i] in ["mRNA", "tRNA", "ncRNA", "rRNA"]:
                if v["strand"] == "+":
                    sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0])
                else:
                    sortedExons = sorted(
                        v["mRNA"][i], key=lambda tup: tup[0], reverse=True
                    )
                annotation[k]["mRNA"][i] = sortedExons
                mrnaSeq = getSeqRegions(SeqRecords, v["contig"], sortedExons)
                annotation[k]["transcript"].append(mrnaSeq)
            if v["type"][i] == "mRNA":
                if v["strand"] == "+":
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0])
                else:
                    sortedCDS = sorted(
                        v["CDS"][i], key=lambda tup: tup[0], reverse=True
                    )
                cdsSeq = getSeqRegions(SeqRecords, v["contig"], sortedCDS)
                # If the translation starts in middle of a codon,
                # we need to truncate the CDS seq either at start or end
                # depending on strand.
                if v["codon_start"][i] > 1:
                    if v["strand"] == "+":
                        # drop first N bases based on codon_start
                        # to reflect the translation frame
                        cdsSeq = cdsSeq[v["codon_start"][i] - 1 :]
                    elif v["strand"] == "-":
                        # drop last N bases based on codon_start
                        # to reflect the translation frame (because this is
                        # is reverse strand gene)
                        endTrunc = len(cdsSeq) - (v["codon_start"][i] - 1)
                        cdsSeq = cdsSeq[0:endTrunc]
                    else:
                        # could trigger more of a warning/error
                        errors.append(
                            "ERROR strand (%s) is nonsensical for %s" % (v["strand"], k)
                        )
                annotation[k]["cds_transcript"].append(cdsSeq)
                annotation[k]["CDS"][i] = sortedCDS
                protSeq = translate(cdsSeq, v["strand"], 0, table=table)
                if protSeq:
                    annotation[k]["protein"].append(protSeq)
                    if protSeq.endswith("*"):
                        annotation[k]["partialStop"][i] = False
                    else:
                        annotation[k]["partialStop"][i] = True
                    if v["codon_start"][i] == 1 and protSeq.startswith("M"):
                        annotation[k]["partialStart"][i] = False
                    else:
                        annotation[k]["partialStart"][i] = True
                # get UTRs
                try:
                    FiveUTR, ThreeUTR = findUTRs(sortedCDS, sortedExons, v["strand"])
                    annotation[k]["5UTR"].append(FiveUTR)
                    annotation[k]["3UTR"].append(ThreeUTR)
                except ValueError:
                    errors.append("ERROR", k, v)
    if debug:
        if len(errors) > 0:
            for e in errors:
                sys.stderr.write("{}\n".format(e))
    return annotation, errors


def findUTRs(cds, mrna, strand):
    FiveUTR = []
    ThreeUTR = []
    if cds != mrna:
        inter = InterLap()
        inter.add(cds)
        for i, x in enumerate(mrna):
            if not x in inter:
                loc = (list(inter)[0][0], list(inter)[-1][1])
                diff = np.subtract(x, loc)
                if diff[0] < 0 and diff[1] < 0:
                    if strand == "+":
                        FiveUTR.append(x)
                    else:
                        ThreeUTR.append(x)
                elif diff[0] > 0 and diff[1] > 0:
                    if strand == "+":
                        ThreeUTR.append(x)
                    else:
                        FiveUTR.append(x)
            else:
                hit = list(inter.find(x))
                if x == hit[0]:
                    continue
                else:
                    diff = np.subtract(x, hit[0])
                    if strand == "+":
                        if int(diff[0]) < 1 and int(diff[1]) == 0:
                            FiveUTR.append((x[0], hit[0][0] - 1))
                        elif int(diff[1]) > 1 and int(diff[0]) == 0:
                            ThreeUTR.append((hit[0][1] + 1, x[1]))
                        elif int(diff[0]) < 1 and int(diff[1]) > 1:
                            FiveUTR.append((x[0], hit[0][0] - 1))
                            ThreeUTR.append((hit[0][1] + 1, x[1]))
                    else:
                        if diff[0] == 0 and diff[1] > 0:
                            FiveUTR.append((hit[0][1] + 1, x[1]))
                        elif diff[0] < 0 and diff[1] == 0:
                            ThreeUTR.append((x[0], hit[0][0] - 1))
                        elif diff[0] < 0 and diff[1] > 0:
                            FiveUTR.append((hit[0][1] + 1, x[1]))
                            ThreeUTR.append((x[0], hit[0][0] - 1))
    return FiveUTR, ThreeUTR


def dict2tbl(
    genesDict,
    scaff2genes,
    scaffLen,
    SeqCenter,
    SeqRefNum,
    skipList,
    output=False,
    annotations=False,
    external=False,
):
    """
    function to take funannotate annotation dictionaries and convert to NCBI tbl output
    """
    duplicates = 0
    pseudo = 0
    nocds = 0
    errors = []
    # to parse annotations, will need to have access to GO OBO dictionary
    goDict = {}
    if annotations:
        goDict, go_format, go_date = go_term_dict()

    def _goFormat(id, goDict=goDict):
        # go_function    serine-type endopeptidase activity|0004252||IEA
        # go_process proteolysis|0006508||IEA
        # go_component   nucleus|0005634||IEA
        if id in goDict:
            if goDict[id]["namespace"] == "biological_process":
                base = "go_process"
            elif goDict[id]["namespace"] == "molecular_function":
                base = "go_function"
            elif goDict[id]["namespace"] == "cellular_component":
                base = "go_component"
            reformatted = "\t\t\t{:}\t{:}|{:}||IEA".format(
                base, goDict[id]["name"], id.replace("GO:", "")
            )
            return reformatted
        else:
            return False

    if output:
        tbl = zopen(output, mode="w")
    else:
        tbl = sys.stdout
    # now convert to tbl format
    for k, v in natsorted(list(scaff2genes.items())):
        tbl.write(">Feature %s\n" % k)
        tbl.write("1\t%s\tREFERENCE\n" % scaffLen.get(k))
        tbl.write("\t\t\t%s\t%s\n" % (SeqCenter, SeqRefNum))
        for genes in v:  # now loop through each gene on the scaffold
            if genes in skipList:
                continue
            # single funannotate standard dictionary
            geneInfo = genesDict.get(genes)
            if not geneInfo["ids"]:
                continue
            if (
                not len(geneInfo["ids"])
                == len(geneInfo["mRNA"])
                == len(geneInfo["CDS"])
            ):
                continue
            #print(genes, geneInfo['note'])
            # check for partial models
            if True in geneInfo["partialStart"]:
                ps = "<"
            else:
                ps = ""
            if True in geneInfo["partialStop"]:
                pss = ">"
            else:
                pss = ""
            # now write gene model
            if geneInfo["strand"] == "+":
                tbl.write(
                    "%s%i\t%s%i\tgene\n"
                    % (ps, geneInfo["location"][0], pss, geneInfo["location"][1])
                )
                if annotations:
                    if geneInfo["name"]:
                        tbl.write("\t\t\tgene\t%s\n" % geneInfo["name"])
                    if geneInfo["gene_synonym"]:
                        for alias in geneInfo["gene_synonym"]:
                            tbl.write("\t\t\tgene_synonym\t%s\n" % alias)
                tbl.write("\t\t\tlocus_tag\t%s\n" % genes)
            else:
                tbl.write(
                    "%s%i\t%s%i\tgene\n"
                    % (ps, geneInfo["location"][1], pss, geneInfo["location"][0])
                )
                if annotations:
                    if geneInfo["name"]:
                        tbl.write("\t\t\tgene\t%s\n" % geneInfo["name"])
                    if geneInfo["gene_synonym"]:
                        for alias in geneInfo["gene_synonym"]:
                            tbl.write("\t\t\tgene_synonym\t%s\n" % alias)
                tbl.write("\t\t\tlocus_tag\t%s\n" % genes)

            # now will output the gene models with -T1, -T2, -T3 annotations based on expression values
            # means need to get the order
            order = []
            # multiple transcripts, so get order of highest TPM
            if len(geneInfo["ids"]) > 1:
                tpms = []
                for num, tpm in enumerate(geneInfo["note"]):
                    for item in tpm:
                        if item.startswith("TPM:"):
                            value = float(item.split(":")[-1])
                            tpms.append((value, num))
                if len(tpms) > 0:
                    for x in sorted(tpms, reverse=True):
                        order.append(x[1])
                else:
                    order = list(range(0, len(geneInfo["ids"])))
            else:
                order.append(0)
            for num, i in enumerate(order):  # now write mRNA and CDS features
                if external:
                    protein_id = geneInfo["ids"][i]
                else:
                    protein_id = genes + "-T" + str(num + 1)
                if geneInfo["type"][i] in ["mRNA", "ncRNA"]:
                    ps = "<"
                    pss = ">"
                    if geneInfo["partialStart"][i] is False:
                        ps = ""
                    if geneInfo["partialStop"][i] is False:
                        pss = ""
                    if geneInfo["strand"] == "+":
                        for num, exon in enumerate(geneInfo["mRNA"][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo["mRNA"][i]) - 1:
                                tbl.write(
                                    "%s%s\t%s%s\t%s\n"
                                    % (ps, exon[0], pss, exon[1], geneInfo["type"][i])
                                )
                            elif num == 0:
                                tbl.write(
                                    "%s%s\t%s\t%s\n"
                                    % (ps, exon[0], exon[1], geneInfo["type"][i])
                                )
                            # this is last one
                            elif num == len(geneInfo["mRNA"][i]) - 1:
                                tbl.write("%s\t%s%s\n" % (exon[0], pss, exon[1]))
                            else:
                                tbl.write("%s\t%s\n" % (exon[0], exon[1]))
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                        tbl.write(
                            "\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n" % (protein_id)
                        )
                        if geneInfo["type"][i] == "mRNA":
                            tbl.write("\t\t\tprotein_id\tgnl|ncbi|%s\n" % (protein_id))
                            for num, cds in enumerate(geneInfo["CDS"][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo["CDS"][i]) - 1:
                                    tbl.write(
                                        "%s%s\t%s%s\tCDS\n" % (ps, cds[0], pss, cds[1])
                                    )
                                elif num == 0:
                                    tbl.write("%s%s\t%s\tCDS\n" % (ps, cds[0], cds[1]))
                                # this is last one
                                elif num == len(geneInfo["CDS"][i]) - 1:
                                    tbl.write("%s\t%s%s\n" % (cds[0], pss, cds[1]))
                                else:
                                    tbl.write("%s\t%s\n" % (cds[0], cds[1]))
                            tbl.write(
                                "\t\t\tcodon_start\t%i\n" % geneInfo["codon_start"][i]
                            )
                            if annotations:  # write functional annotation
                                if geneInfo["EC_number"][i]:
                                    for EC in geneInfo["EC_number"][i]:
                                        tbl.write("\t\t\tEC_number\t%s\n" % EC)
                                if geneInfo["db_xref"][i]:
                                    for xref in geneInfo["db_xref"][i]:
                                        tbl.write("\t\t\tdb_xref\t%s\n" % xref)
                                if geneInfo["go_terms"][i]:
                                    for go in geneInfo["go_terms"][i]:
                                        goLine = _goFormat(go)
                                        if goLine:
                                            tbl.write("{:}\n".format(goLine))
                                if geneInfo["note"][i]:
                                    for item in geneInfo["note"][i]:
                                        tbl.write("\t\t\tnote\t%s\n" % item)
                            tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                            tbl.write(
                                "\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n" % (protein_id)
                            )
                            tbl.write("\t\t\tprotein_id\tgnl|ncbi|%s\n" % (protein_id))
                    else:  # means this is on crick strand
                        for num, exon in enumerate(geneInfo["mRNA"][i]):
                            # single exon, so slightly differnt method
                            if num == 0 and num == len(geneInfo["mRNA"][i]) - 1:
                                tbl.write(
                                    "%s%s\t%s%s\t%s\n"
                                    % (ps, exon[1], pss, exon[0], geneInfo["type"][i])
                                )
                            elif num == 0:
                                tbl.write(
                                    "%s%s\t%s\t%s\n"
                                    % (ps, exon[1], exon[0], geneInfo["type"][i])
                                )
                            # this is last one
                            elif num == len(geneInfo["mRNA"][i]) - 1:
                                tbl.write("%s\t%s%s\n" % (exon[1], pss, exon[0]))
                            else:
                                tbl.write("%s\t%s\n" % (exon[1], exon[0]))
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                        tbl.write(
                            "\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n" % (protein_id)
                        )
                        if geneInfo["type"][i] == "mRNA":
                            tbl.write("\t\t\tprotein_id\tgnl|ncbi|%s\n" % (protein_id))
                            for num, cds in enumerate(geneInfo["CDS"][i]):
                                # single exon, so slightly differnt method
                                if num == 0 and num == len(geneInfo["CDS"][i]) - 1:
                                    tbl.write(
                                        "%s%s\t%s%s\tCDS\n" % (ps, cds[1], pss, cds[0])
                                    )
                                elif num == 0:
                                    tbl.write("%s%s\t%s\tCDS\n" % (ps, cds[1], cds[0]))
                                # this is last one
                                elif num == (len(geneInfo["CDS"][i]) - 1):
                                    tbl.write("%s\t%s%s\n" % (cds[1], pss, cds[0]))
                                else:
                                    tbl.write("%s\t%s\n" % (cds[1], cds[0]))
                            tbl.write(
                                "\t\t\tcodon_start\t%i\n" % geneInfo["codon_start"][i]
                            )
                            if annotations:  # write functional annotation
                                if geneInfo["EC_number"][i]:
                                    for EC in geneInfo["EC_number"][i]:
                                        tbl.write("\t\t\tEC_number\t%s\n" % EC)
                                if geneInfo["db_xref"][i]:
                                    for xref in geneInfo["db_xref"][i]:
                                        tbl.write("\t\t\tdb_xref\t%s\n" % xref)
                                if geneInfo["go_terms"][i]:
                                    for go in geneInfo["go_terms"][i]:
                                        goLine = _goFormat(go)
                                        if goLine:
                                            tbl.write("{:}\n".format(goLine))
                                if geneInfo["note"][i]:
                                    for item in geneInfo["note"][i]:
                                        tbl.write("\t\t\tnote\t%s\n" % item)
                            tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                            tbl.write(
                                "\t\t\ttranscript_id\tgnl|ncbi|%s_mrna\n" % (protein_id)
                            )
                            tbl.write("\t\t\tprotein_id\tgnl|ncbi|%s\n" % (protein_id))
                elif geneInfo["type"][i] == "tRNA":
                    if geneInfo["strand"] == "+":
                        for num, exon in enumerate(geneInfo["mRNA"][i]):
                            if num == 0:
                                tbl.write(
                                    "%s\t%s\t%s\n"
                                    % (exon[0], exon[1], geneInfo["type"][i])
                                )
                            else:
                                tbl.write("%s\t%s\n" % (exon[0], exon[1]))
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                        if geneInfo["product"] == "tRNA-Xxx":
                            tbl.write("\t\t\tpseudo\n")
                    else:
                        for num, exon in enumerate(geneInfo["mRNA"][i]):
                            if num == 0:
                                tbl.write(
                                    "%s\t%s\t%s\n"
                                    % (exon[1], exon[0], geneInfo["type"][i])
                                )
                            else:
                                tbl.write("%s\t%s\n" % (exon[1], exon[0]))
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                        if geneInfo["product"] == "tRNA-Xxx":
                            tbl.write("\t\t\tpseudo\n")
                elif geneInfo["type"][i] in ["rRNA"]:
                    if geneInfo["strand"] == "+":
                        tbl.write(
                            "%s\t%s\t%s\n"
                            % (
                                geneInfo["location"][0],
                                geneInfo["location"][1],
                                geneInfo["type"][i],
                            )
                        )
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
                    else:
                        tbl.write(
                            "%s\t%s\t%s\n"
                            % (
                                geneInfo["location"][1],
                                geneInfo["location"][0],
                                geneInfo["type"][i],
                            )
                        )
                        tbl.write("\t\t\tproduct\t%s\n" % geneInfo["product"][i])
    if output:
        tbl.close()
    return errors, duplicates, pseudo, nocds

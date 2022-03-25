import sys
from collections import defaultdict, OrderedDict
from natsort import natsorted
from itertools import product
import numpy as np
from .utils import zopen
from .gff import gff2dict
from .interlap import InterLap
from .consensus import getAED


def compare(args):
    compareAnnotations(args.reference, args.query, args.fasta, output=args.out)


def gff2interlap(input, fasta):
    """
    function to parse GFF3 file, construct scaffold/gene interlap dictionary and funannotate standard annotation dictionary
    """
    inter = defaultdict(InterLap)
    Genes = gff2dict(input, fasta)
    for k, v in natsorted(list(Genes.items())):
        inter[v["contig"]].add((v["location"][0], v["location"][1], k))
    return inter, Genes


def countFeatures(input):
    # given funannotate dictionary, count up some general features
    mRNAg, mRNAt, tRNAg, tRNAt = (0,) * 4
    for k, v in natsorted(list(input.items())):
        if v["type"] == "mRNA":
            mRNAg += 1
            mRNAt += len(v["ids"])
        elif v["type"] == "tRNA":
            tRNAg += 1
            tRNAt += len(v["ids"])
    return len(input), mRNAg, mRNAt, tRNAg, tRNAt


def pairwiseAED(query, reference):
    """
    takes a multiple transcripts and sums AED from lowest pairwise comparison and then calculates
    the average based on number of transcripts in the query
    """
    AEDsum = []
    pAED = [float(getAED(a, b)) for a, b in product(query, reference)]
    # split into parts to get lowest AED
    splitAED = [pAED[i : i + len(query)] for i in range(0, len(pAED), len(query))]
    for pair in splitAED:
        AEDsum.append(min(pair))
    AEDavg = sum(AEDsum) / len(query)
    return "{:.3f}".format(AEDavg)


def compareAnnotations(old, new, fasta, output=False):
    """
    function takes two GFF annotated genomes and compares gene models
    output is a tsv file for each locus and a description of what is different
    can handle multiple transcripts per locus
    """
    if output:
        out = zopen(output, mode="w")
    else:
        out = sys.stdout
    result = {}
    global no_change, identicalCDS, refUnique, queryUnique
    no_change, identicalCDS, refUnique, queryUnique, totalmatches, totallength = (
        0,
    ) * 6
    oldInter, oldGenes = gff2interlap(old, fasta)
    NumOldLoci, NumOldGenes, NumOldmRNA, NumOldtRNALoci, NumOldtRNA = countFeatures(
        oldGenes
    )
    sys.stderr.write(
        "Reference: {} contigs {} genes\n".format(len(oldInter), len(oldGenes))
    )
    sys.stderr.write(
        "{}\n".format([NumOldLoci, NumOldGenes, NumOldmRNA, NumOldtRNALoci, NumOldtRNA])
    )
    newInter, newGenes = gff2interlap(new, fasta)
    sys.stderr.write(
        "Query: {} contigs {} genes\n".format(len(newInter), len(newGenes))
    )
    NumNewLoci, NumNewGenes, NumNewmRNA, NumNewtRNALoci, NumNewtRNA = countFeatures(
        newGenes
    )
    sys.stderr.write(
        "{}\n".format([NumNewLoci, NumNewGenes, NumNewmRNA, NumNewtRNALoci, NumNewtRNA])
    )

    # now go through the updated annotation, comparing to old annot
    ref_seen = set()
    for contig in newInter:
        for gene in newInter[contig]:
            # means this is a new model, so add it
            hitList = list(oldInter[contig].find((gene[0], gene[1])))
            if len(hitList) < 1:
                result[gene[2]] = {
                    "contig": newGenes[gene[2]]["contig"],
                    "location": newGenes[gene[2]]["location"],
                    "ref_type": None,
                    "ref_location": None,
                    "query_location": newGenes[gene[2]]["location"],
                    "query_id": gene[2],
                    "query_type": newGenes[gene[2]]["mRNA"],
                    "ref_id": None,
                    "cdsAED": "1.000",
                    "exonAED": "1.000",
                    "ref_transcripts": 0,
                    "query_transcripts": len(newGenes[gene[2]]["ids"]),
                    "ref_strand": None,
                    "query_strand": newGenes[gene[2]]["strand"],
                }
                continue
            # there might be some overlapping transcripts, so get best hit?
            hit = []
            # get best hit
            for z in hitList:
                diffs = np.subtract((gene[0], gene[1]), (z[0], z[1]))
                totaldiffs = abs(diffs[0]) + abs(diffs[1])
                hit.append((totaldiffs, z[2]))
            besthit = min(hit)

            # get the old annotation
            hitInfo = oldGenes.get(besthit[1])
            ref_seen.add(besthit[1])

            # calculate AED
            exonAED = pairwiseAED(newGenes[gene[2]]["mRNA"], hitInfo["mRNA"])
            if "mRNA" in newGenes[gene[2]]["type"] and "mRNA" in hitInfo["type"]:
                cdsAED = pairwiseAED(newGenes[gene[2]]["CDS"], hitInfo["CDS"])
            else:
                cdsAED = "0.000"

            if not besthit[1] in result:
                result[besthit[1]] = {
                    "contig": newGenes[gene[2]]["contig"],
                    "location": hitInfo["location"],
                    "ref_type": hitInfo["type"],
                    "ref_location": hitInfo["location"],
                    "query_location": newGenes[gene[2]]["location"],
                    "query_id": gene[2],
                    "query_type": newGenes[gene[2]]["type"],
                    "cdsAED": cdsAED,
                    "exonAED": exonAED,
                    "ref_transcripts": len(hitInfo["ids"]),
                    "query_transcripts": len(newGenes[gene[2]]["ids"]),
                    "ref_strand": hitInfo["strand"],
                    "query_strand": newGenes[gene[2]]["strand"],
                    "ref_id": besthit[1],
                }
            # get some summary stats as you loop through
            if float(exonAED) == 0 and float(cdsAED) == 0:
                no_change += 1
            elif float(cdsAED) == 0:
                identicalCDS += 1
    # now add old genes that did not have overlaps
    for contig in oldInter:
        for gene in oldInter[contig]:
            if not gene[2] in ref_seen:
                result[gene[2]] = {
                    "contig": oldGenes[gene[2]]["contig"],
                    "location": oldGenes[gene[2]]["location"],
                    "ref_type": oldGenes[gene[2]]["type"][0],
                    "ref_location": oldGenes[gene[2]]["location"],
                    "query_location": None,
                    "query_id": None,
                    "query_type": None,
                    "ref_id": gene[2],
                    "cdsAED": "1.000",
                    "exonAED": "1.000",
                    "ref_transcripts": len(oldGenes[gene[2]]["ids"]),
                    "query_transcripts": 0,
                    "ref_strand": oldGenes[gene[2]]["strand"],
                    "query_strand": None,
                }
    total_cdsAED = []
    total_exonAED = []

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(result.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    out.write(
        "Reference_Location\tReference_ID\tRef_strand\tRef_Num_Transcripts\tQuery_Location\tQuery_ID\tQuery_strand\tQuery_Num_Transcripts\tmRNA_AED\tCDS_AED\n"
    )
    for k, v in list(sortedGenes.items()):
        Rstart = str(v["location"][0])
        Rend = str(v["location"][1])
        if v["query_id"]:
            Qstart = str(v["query_location"][0])
            Qend = str(v["query_location"][1])
        else:
            Qstart = "None"
            Qend = "None"
        total_cdsAED.append(float(v["cdsAED"]))
        total_exonAED.append(float(v["exonAED"]))
        out.write(
            "{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}:{:}-{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                v["contig"],
                Rstart,
                Rend,
                v["ref_id"],
                v["ref_strand"],
                v["ref_transcripts"],
                v["contig"],
                Qstart,
                Qend,
                v["query_id"],
                v["query_strand"],
                v["query_transcripts"],
                v["exonAED"],
                v["cdsAED"],
            )
        )
    Avg_cdsAED = sum(total_cdsAED) / float(len(total_cdsAED))
    Avg_exonAED = sum(total_exonAED) / float(len(total_exonAED))
    totalPident = 0.00
    if totalmatches > 0:
        totalPident = totalmatches / totallength
    if output:
        out.close()
    return [
        NumOldLoci,
        NumOldGenes,
        NumOldmRNA,
        NumOldtRNALoci,
        NumOldtRNA,
        refUnique,
        no_change,
        identicalCDS,
        0.000,
        0.000,
        1,
        NumNewLoci,
        NumNewGenes,
        NumNewmRNA,
        NumNewtRNALoci,
        NumNewtRNA,
        queryUnique,
        no_change,
        identicalCDS,
        Avg_exonAED,
        Avg_cdsAED,
        totalPident,
    ]

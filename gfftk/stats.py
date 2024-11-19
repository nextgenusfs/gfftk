import sys
import json
from statistics import median
from .utils import check_inputs, zopen
from .gff import gff2dict
from .genbank import tbl2dict


def stats(args):
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
            raise SystemExit(1)
    # okay now we can load and convert
    if args.input_format == "tbl":
        Genes, parse_errors = tbl2dict(args.input, args.fasta, {}, table=1)
    elif args.input_format == "gff3":
        Genes = gff2dict(args.input, args.fasta, table=1, debug=args.debug)
    # calculate annotation stats and output
    annot_stats = annotation_stats(Genes)
    if args.out:
        with zopen(args.out, mode="w") as outfile:
            outfile.write(json.dumps(annot_stats, indent=2))
    else:
        sys.stdout.write(json.dumps(annot_stats, indent=2))


def annotation_stats(Genes):
    """
    function to output annotation stats from annotation dictionary
    """
    stats = {
        "genes": 0,
        "common_name": 0,
        "mRNA": 0,
        "tRNA": 0,
        "ncRNA": 0,
        "rRNA": 0,
        "repeat_region": 0,
        "misc_feature": 0,
        "avg_gene_length": 0,
        "transcript-level": {
            "CDS_transcripts": 0,
            "CDS_five_utr": 0,
            "CDS_three_utr": 0,
            "CDS_no_utr": 0,
            "CDS_five_three_utr": 0,
            "CDS_complete": 0,
            "CDS_no-start": 0,
            "CDS_no-stop": 0,
            "CDS_no-start_no-stop": 0,
            "total_exons": 0,
            "total_cds_exons": 0,
            "average_number_transcripts_per_gene": 0,
            "multiple_exon_transcript": 0,
            "single_exon_transcript": 0,
            "average_number_cds_exons": 0,
            "avg_exon_length": 0,
            "median_number_exons": 0,
            "max_number_exons": 0,
            "avg_protein_length": 0,
            "avg_transcript_length": 0,
            "functional": {
                "go_terms": 0,
                "interproscan": 0,
                "eggnog": 0,
                "pfam": 0,
                "cazyme": 0,
                "merops": 0,
                "busco": 0,
                "secretion": 0,
            },
        },
    }
    # parse annotation dictionary
    if len(Genes) > 0:
        protLengths = []
        geneLengths = []
        exonLengths = []
        transcriptLengths = []
        num_transcripts_per_gene = []
        cdsExons = []
        for k, v in Genes.items():
            stats["genes"] += 1
            gLength = v["location"][1] - v["location"][0]
            geneLengths.append(gLength)
            if "tRNA" in v["type"]:
                stats["tRNA"] += 1
            elif "rRNA" in v["type"]:
                stats["rRNA"] += 1
            elif "repeat_region" in v["type"]:
                stats["repeat_region"] += 1
            elif v["type"] == "misc_feature":
                stats["misc_feature"] += 1
            elif "mRNA" not in v["type"] and "ncRNA" in v["type"]:
                stats["ncRNA"] += 1
            if v["name"]:
                stats["common_name"] += 1
            num_transcripts_per_gene.append(len(v["ids"]))
            for i in range(0, len(v["ids"])):
                if v["type"][i] == "mRNA":
                    stats["mRNA"] += 1
                    stats["transcript-level"]["CDS_transcripts"] += 1
                    pLen = len(v["protein"][i])
                    if v["protein"][i].endswith("*"):
                        pLen -= 1
                    protLengths.append(pLen)
                    transcriptLengths.append(sum([x[1] - x[0] for x in v["mRNA"][i]]))
                    if len(v["mRNA"][i]) > 1:
                        stats["transcript-level"]["multiple_exon_transcript"] += 1
                    else:
                        stats["transcript-level"]["single_exon_transcript"] += 1
                    try:
                        for y in v["mRNA"][i]:
                            exon_length = y[1] - y[0]
                            exonLengths.append(exon_length)
                    except IndexError:
                        sys.stderr.write(
                            "ERROR calculating exon length for {}:\n{}".format(
                                k, json.dumps(v, indent=2)
                            )
                        )
                        raise SystemExit(1)
                    stats["transcript-level"]["total_exons"] += len(v["mRNA"][i])
                    stats["transcript-level"]["total_exons"] += len(v["5UTR"][i])
                    stats["transcript-level"]["total_exons"] += len(v["3UTR"][i])
                    stats["transcript-level"]["total_cds_exons"] += len(v["CDS"][i])
                    cdsExons.append(len(v["CDS"][i]))
                    if v["partialStart"][i] and v["partialStop"][i]:
                        stats["transcript-level"]["CDS_no-start_no-stop"] += 1
                    elif v["partialStart"][i]:
                        stats["transcript-level"]["CDS_no-start"] += 1
                    elif v["partialStop"][i]:
                        stats["transcript-level"]["CDS_no-stop"] += 1
                    else:
                        stats["transcript-level"]["CDS_complete"] += 1
                    if len(v["5UTR"][i]) > 0 and len(v["3UTR"][i]) > 0:
                        stats["transcript-level"]["CDS_five_three_utr"] += 1
                    elif len(v["3UTR"][i]) > 0:
                        stats["transcript-level"]["CDS_three_utr"] += 1
                    elif len(v["5UTR"][i]) > 0:
                        stats["transcript-level"]["CDS_three_utr"] += 1
                    else:
                        stats["transcript-level"]["CDS_no_utr"] += 1
                    if v["go_terms"][i]:
                        stats["transcript-level"]["functional"]["go_terms"] += 1
                    if any(s.startswith("PFAM:") for s in v["db_xref"][i]):
                        stats["transcript-level"]["functional"]["pfam"] += 1
                    if any(s.startswith("InterPro:") for s in v["db_xref"][i]):
                        stats["transcript-level"]["functional"]["interproscan"] += 1
                    if any(s.startswith("EggNog:") for s in v["note"][i]):
                        stats["transcript-level"]["functional"]["eggnog"] += 1
                    if any(s.startswith("CAZy:") for s in v["note"][i]):
                        stats["transcript-level"]["functional"]["cazyme"] += 1
                    if any(s.startswith("MEROPS:") for s in v["note"][i]):
                        stats["transcript-level"]["functional"]["merops"] += 1
                    if any(s.startswith("BUSCO:") for s in v["note"][i]):
                        stats["transcript-level"]["functional"]["busco"] += 1
                    if any(s.startswith("SECRETED:") for s in v["note"][i]):
                        stats["transcript-level"]["functional"]["secretion"] += 1

        stats["avg_gene_length"] = round(sum(geneLengths) / float(len(geneLengths)), 2)
        stats["transcript-level"]["avg_protein_length"] = round(
            sum(protLengths) / float(len(protLengths)), 2
        )
        try:
            stats["transcript-level"]["avg_exon_length"] = round(
                sum(exonLengths) / float(len(exonLengths)), 2
            )
        except ZeroDivisionError:
            stats["transcript-level"]["avg_exon_length"] = 0
        stats["transcript-level"]["average_number_transcripts_per_gene"] = round(
            sum(num_transcripts_per_gene) / float(len(num_transcripts_per_gene)), 2
        )
        stats["transcript-level"]["avg_transcript_length"] = round(
            sum(transcriptLengths) / float(len(transcriptLengths)), 2
        )
        stats["transcript-level"]["average_number_cds_exons"] = round(
            sum(cdsExons) / float(len(cdsExons)), 2
        )
        stats["transcript-level"]["max_number_exons"] = max(cdsExons)
        stats["transcript-level"]["median_number_exons"] = median(cdsExons)
    return stats

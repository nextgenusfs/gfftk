import itertools
import json
import os
import random
import sys
from collections import Counter, OrderedDict, defaultdict
import numpy as np
from natsort import natsorted
import uuid
from .fasta import fastaparser
from .gff import gff2dict
from . import interlap
from .log import startLogging, system_info
from .utils import check_inputs, zopen


def consensus(args):
    logger = startLogging(logfile=args.logfile)
    global log
    if args.silent:
        log = logger.debug
    else:
        log = logger.info
    system_info(log)
    logger.info(args)
    check_inputs(
        [args.fasta] + args.genes + args.proteins + args.transcripts + args.repeats
    )
    _ = generate_consensus(
        args.fasta,
        args.genes,
        args.proteins,
        args.transcripts,
        args.weights,
        args.out,
        debug=args.debug,
        minscore=args.minscore,
        repeats=args.repeats,
        repeat_overlap=args.repeat_overlap,
        min_exon=args.min_exon,
        max_exon=args.max_exon,
        min_intron=args.min_intron,
        max_intron=args.max_intron,
        log=log,
    )


def generate_consensus(
    fasta,
    genes,
    proteins,
    transcripts,
    weights,
    out,
    debug=False,
    minscore=False,
    repeats=False,
    repeat_overlap=90,
    tiebreakers="calculated",
    min_exon=3,
    min_intron=10,
    max_intron=-1,
    max_exon=-1,
    evidence_derived_models=[],
    log=sys.stderr.write,
):
    log("GFFtk consensus will generate the best gene model at each locus")
    if debug:
        debug = out + ".all_gene_models.gff3"
        if os.path.isfile(debug):
            os.remove(debug)
    if weights:
        WEIGHTS = {}
        for x in weights:
            if ":" in x:
                source, w = x.split(":", 1)
            else:
                source = x
                w = 1
            WEIGHTS[source] = int(w)
    else:
        WEIGHTS = {}
    log("Parsing GFF3 files and clustering data into strand specific loci")
    data = parse_data(fasta, genes, proteins, transcripts, log=log)
    if tiebreakers == "calculated":
        # this will use evidence to calculate source tiebreakers
        order, n_evidence = calculate_source_order(data)
        if n_evidence > 0:
            log(
                "Filtered gene models for evidence: {} loci have >3 genes and >2 alignment evidence".format(
                    n_evidence
                )
            )
    else:
        tmpD = {}
        for contig, obj in data.items():
            for strand in ["+", "-"]:
                for locus in obj[strand]:
                    for gene in locus["genes"]:
                        name, source, coords, cstart = gene
                        if source in WEIGHTS:
                            ord_score = WEIGHTS.get(source)
                        else:
                            ord_score = 1
                        if source not in tmpD:
                            tmpD[source] = ord_score
        order = OrderedDict(sorted(tmpD.items(), key=lambda kv: kv[1], reverse=True))
    log(
        "Using these filtered loci, the calculated gene model source weights to use as tiebreakers:\n{}".format(
            json.dumps(order)
        )
    )

    consensus = {}
    plus_consensus = defaultdict(interlap.InterLap)
    minus_consensus = defaultdict(interlap.InterLap)
    counter = 1
    if debug:
        locus_bed = zopen(out + ".loci.bed", mode="w")
    else:
        locus_bed = open(os.devnull, "w")
    split_stats = {}
    for contig, obj in data.items():
        for strand in ["+", "-"]:
            for locus in obj[strand]:
                name = "locus_{}".format(counter)
                # print(name, contig, strand, locus["locus"])
                # keepers is a list of lists
                keepers = best_model(
                    name,
                    contig,
                    strand,
                    locus,
                    weights=WEIGHTS,
                    order=order,
                    debug=debug,
                    min_exon=min_exon,
                    min_intron=min_intron,
                    max_intron=max_intron,
                    max_exon=max_exon,
                    evidence_derived_models=evidence_derived_models,
                )
                # print(name, keepers)
                if not keepers:
                    continue
                if len(keepers) not in split_stats:
                    split_stats[len(keepers)] = 1
                else:
                    split_stats[len(keepers)] += 1
                # write locus bed file
                locus_bed.write(
                    f'{contig}\t{locus["locus"][0]}\t{locus["locus"][1]}\t{name};n_genes={len(keepers)}\t0\t{strand}\n'
                )
                if contig not in consensus:
                    consensus[contig] = {}
                for keep in keepers:
                    if keep[0] not in consensus[contig]:
                        consensus[contig][keep[0]] = {
                            "strand": strand,
                            "locus": name,
                            "codon_start": keep[1]["codon_start"],
                            "coords": keep[1]["coords"],
                            "score": keep[1]["score"],
                            "source": keep[1]["source"],
                        }
                        if strand == "+":
                            plus_consensus[contig].add(
                                (
                                    min(min(keep[1]["coords"])),
                                    max(max(keep[1]["coords"])),
                                    name,
                                    keep[1]["score"],
                                )
                            )
                        elif strand == "-":
                            minus_consensus[contig].add(
                                (
                                    min(min(keep[1]["coords"])),
                                    max(max(keep[1]["coords"])),
                                    name,
                                    keep[1]["score"],
                                )
                            )
                counter += 1
    locus_bed.close()
    # debug models here
    # log(split_stats)
    # check for minimum score
    if not minscore:
        score_threshold = auto_score_threshold(WEIGHTS, order)
    else:
        score_threshold = int(minscore)
    log("Setting minimum gene model score to {}".format(score_threshold))
    sources = {}
    filtered = {}
    total = 0
    for c, o in consensus.items():
        for m, v in o.items():
            if v["score"] > score_threshold:
                # now see if model is contained on opposite strand, if so then skip
                if v["strand"] == "+":
                    hits = list(
                        minus_consensus[c].find(
                            ((min(min(v["coords"])), max(max(v["coords"]))))
                        )
                    )
                elif v["strand"] == "-":
                    hits = list(
                        plus_consensus[c].find(
                            ((min(min(v["coords"])), max(max(v["coords"]))))
                        )
                    )
                if len(hits) > 0:
                    contain = [
                        contained(
                            (min(min(v["coords"])), max(max(v["coords"]))), (x[0], x[1])
                        )
                        for x in hits
                    ]
                    if any(contain):
                        continue
                total += 1
                if v["source"] not in sources:
                    sources[v["source"]] = 1
                else:
                    sources[v["source"]] += 1
                v["contig"] = c
                v["name"] = m
                v["location"] = (min(min(v["coords"])), max(max(v["coords"])))
                filtered[m] = v

    # now we can filter models in repeat regions
    if repeats:
        final, dropped_n = filter_models_repeats(
            fasta, repeats, filtered, filter_threshold=repeat_overlap, log=log
        )
        log("{} gene models were dropped due to repeat overlap".format(dropped_n))
    else:
        final = filtered
    log(
        "{} consensus gene models derived from these sources:\n{}".format(
            len(final),
            json.dumps(sorted(sources.items(), key=lambda x: x[1], reverse=True)),
        )
    )
    # finally we can write to GFF3 format
    gff_writer(final, out)
    log("GFFtk consensus is finished: {}".format(out))
    return final


def filter_models_repeats(fasta, repeats, gene_models, filter_threshold=90, log=False):
    dropped = 0
    # return fasta length for stats generation
    seq_length = fasta_length(fasta)
    # build interlap object with repeats
    repeat_length = 0
    repeat_inter = defaultdict(interlap.InterLap)
    for r in repeats:
        if os.path.isfile(r):
            if r.endswith(".bed"):
                repeat_inter, repeat_length = bed2interlap(r, inter=repeat_inter)
            else:
                repeat_inter, repeat_length = gff2interlap(r, fasta, inter=repeat_inter)
    pct_repeats = repeat_length / seq_length
    if log:
        log(
            "Loaded repeats representing {:.2%} of the genome and filtering out loci that are > {}% overlap with repeats".format(
                pct_repeats, filter_threshold
            )
        )
    # now we can filter results by going to through gene models and
    # looking at which models overlap with repeats
    filtered = {}
    for k, v in gene_models.items():
        hits = list(repeat_inter[v["contig"]].find(v["location"]))
        if len(hits) > 0:
            gene_length = v["location"][1] - v["location"][0]
            repeat_overlap = 0
            for x in hits:
                if v["location"][0] > x[0]:
                    s = v["location"][0]
                else:
                    s = x[0]
                if v["location"][1] < x[1]:
                    e = v["location"][1]
                else:
                    e = x[1]
                repeat_overlap += e - s
            gene_repeat_overlap = repeat_overlap / gene_length
            if gene_repeat_overlap <= (filter_threshold / 100):
                filtered[k] = v
            else:
                dropped += 1
        else:
            filtered[k] = v
    return filtered, dropped


def fasta_length(fasta):
    length = 0
    with zopen(fasta) as infile:
        for title, seq in fastaparser(infile):
            length += len(seq)
    return length


def gff2interlap(input, fasta, inter=False):
    """
    function to parse GFF3 file, construct scaffold/gene interlap dictionary
    """
    length = 0
    if not inter:
        inter = defaultdict(interlap.InterLap)
    Genes = gff2dict(input, fasta)
    for k, v in natsorted(list(Genes.items())):
        inter[v["contig"]].add((v["location"][0], v["location"][1]))
        length += v["location"][0] - v["location"][1]
    return inter, length


def bed2interlap(bedfile, inter=False):
    # load interlap object from a bed file
    length = 0
    if not inter:
        inter = defaultdict(interlap.InterLap)
    with zopen(bedfile) as infile:
        for line in infile:
            line = line.strip()
            chr, start, end = line.split("\t")[:3]
            inter[chr].add((int(start), int(end)))
            length += int(end) - int(start)
    return inter, length


def auto_score_threshold(weights, order, user_weight=6):
    # figure out minimum score
    # min_weights = min(weights.values())
    # min_order = min(order.values())
    allweights = {}
    for w in order.keys():
        allweights[w] = weights.get(w, 1) * user_weight
    min_weights = min(allweights.values())
    threshold = 1 + min_weights
    return threshold


def cluster_interlap(obj):
    # use the interlap.reduce function to get clusters
    # then come back through and assign data to the clusters
    all_coords = [(x[0], x[1]) for x in obj]
    bins = interlap.reduce(all_coords)
    results = []
    for b in bins:
        hits = list(obj.find(b))
        if len(hits) > 1:
            results.append(
                {
                    "locus": b,
                    "genes": [x[2:] for x in hits],
                    "proteins": [],
                    "transcripts": [],
                    "repeats": [],
                }
            )
    return results


def cluster_interlap_old(obj):
    seen = set()
    results = []
    for i, loc in enumerate(obj):
        if loc[2] in seen:
            continue
        hits = list(obj.find((loc[0], loc[1])))
        if len(hits) > 1:
            for x in hits:
                seen.add(x[2])
            # sub cluster for unique sources
            cleaned = [hits]
            for r in cleaned:
                starts = [x[0] for x in r]
                ends = [x[1] for x in r]
                data = [x[2:] for x in r]
                results.append(
                    {
                        "locus": (min(starts), max(ends)),
                        "genes": data,
                        "proteins": [],
                        "transcripts": [],
                        "repeats": [],
                    }
                )
        else:
            results.append(
                {
                    "locus": (loc[0], loc[1]),
                    "genes": [loc[2:]],
                    "proteins": [],
                    "transcripts": [],
                    "repeats": [],
                }
            )
    return results


def sub_cluster(obj):
    # input will be a list of lists for each locus
    # need to check the source and split into sub clusters if any source > 1
    sources = [x[3] for x in obj]
    s, n = Counter(sources).most_common(1)[0]
    if n == 1:  # then no duplicate sources, return
        return [obj]
    else:  # now we need to split these up
        seeds = []
        queries = []
        for x in obj:
            if x[3] == s:
                seeds.append([x])
            else:
                queries.append(x)
        # now we can assign the queries to seeds
        for y in queries:
            scores = []
            for w in seeds:
                o = get_overlap([y[0], y[1]], [w[0][0], w[0][1]])
                scores.append(o)
            idx = scores.index(max(scores))
            seeds[idx].append(y)
        return seeds


def contained(a, b):
    # check if coords in a are completely contained in b
    left = a[0] - b[0]
    right = a[1] - b[1]
    if left > 0 and right < 0:
        return True
    else:
        return False


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_loci(annot_dict):
    # input is annotation dictionary
    plus_inter = defaultdict(interlap.InterLap)
    minus_inter = defaultdict(interlap.InterLap)
    contigs = []
    pseudo = []
    for gene, info in annot_dict.items():
        if info["contig"] not in contigs:
            contigs.append(info["contig"])
        if info["pseudo"] is True:
            pseudo.append((gene, info))
            continue
        if info["protein"][0].count("*") > 1:
            pseudo.append((gene, info))
            continue
        if info["strand"] == "+":
            plus_inter[info["contig"]].add(
                [
                    info["location"][0],
                    info["location"][1],
                    gene,
                    info["source"],
                    info["CDS"],
                    info["codon_start"][0],
                ]
            )
        else:
            minus_inter[info["contig"]].add(
                [
                    info["location"][0],
                    info["location"][1],
                    gene,
                    info["source"],
                    info["CDS"],
                    info["codon_start"][0],
                ]
            )
    loci = {}
    n_loci = 0
    for chr in contigs:
        plus_bins = {}
        minus_bins = {}
        if chr in plus_inter:
            plus_bins = cluster_interlap(plus_inter[chr])
        n_loci += len(plus_bins)
        if chr in minus_inter:
            minus_bins = cluster_interlap(minus_inter[chr])
        n_loci += len(minus_bins)
        loci[chr] = {"+": plus_bins, "-": minus_bins}
    return loci, n_loci, pseudo


def gffevidence2dict(file, Evi):
    """
    general function to take a evidence in GFF3 format and return locustag/ID dictionary
    """
    seen = set()
    with zopen(file, mode="r") as input:
        for line in input:
            if line.startswith("\n") or line.startswith("#"):
                continue
            line = line.rstrip()
            if line in seen:
                continue
            try:
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
            except:  # noqa: E722
                continue
            seen.add(line)
            start = int(start)
            end = int(end)
            info = attributes.split(";")
            fields = {}
            for field in info:
                try:
                    k, v = field.split("=")
                    fields[k] = v
                except (IndexError, ValueError) as E:
                    pass
            if "ID" not in fields:
                continue
            ID = fields["ID"]
            if "target" in fields:
                Target = fields["target"]
            else:
                Target = None
            if ID not in Evi:
                Evi[ID] = {
                    "target": Target,
                    "type": feature,
                    "source": source,
                    "strand": strand,
                    "phase": [phase],
                    "contig": contig,
                    "coords": [(start, end)],
                    "location": (start, end),
                    "score": [score],
                }
            else:
                Evi[ID]["coords"].append((start, end))
                Evi[ID]["score"].append(score)
                Evi[ID]["phase"].append(phase)
                if start < Evi[ID]["location"][0]:
                    Evi[ID]["location"] = (start, Evi[ID]["location"][1])
                if end > Evi[ID]["location"][1]:
                    Evi[ID]["location"] = (Evi[ID]["location"][0], end)
    return Evi


def add_evidence(loci, evidence, source="proteins"):
    # here we will build interlap obj of evidence
    # then loop through loci and pull in evidence that aligns to each locus
    plus_inter = defaultdict(interlap.InterLap)
    minus_inter = defaultdict(interlap.InterLap)
    for k, v in evidence.items():
        if v["strand"] == "+":
            plus_inter[v["contig"]].add(
                [v["location"][0], v["location"][1], k, v["source"], v["coords"]]
            )
        else:
            minus_inter[v["contig"]].add(
                [v["location"][0], v["location"][1], k, v["source"], v["coords"]]
            )
    for k, v in loci.items():
        for i, x in enumerate(v["+"]):
            hits = list(plus_inter[k].find(x["locus"]))
            if len(hits) > 0:
                # print(k, '+', x['locus'], i, hits)
                cleaned = []
                for h in hits:
                    cleaned.append(h[2:])
                loci[k]["+"][i][source] = cleaned
        for z, y in enumerate(v["-"]):
            hits2 = list(minus_inter[k].find(y["locus"]))
            if len(hits) > 0:
                # print(k, '-', x['locus'], z, hits)
                cleaned = []
                for h in hits2:
                    cleaned.append(h[2:])
                loci[k]["-"][z][source] = cleaned


def ensure_unique_names(genes):
    # takes annotation dictionary and appends unique slug
    slug = str(uuid.uuid4())[:8]
    cleaned = {}
    for k, v in genes.items():
        cleaned[f"{k}.{slug}"] = v
    return cleaned


def parse_data(genome, gene, protein, transcript, log=sys.stderr.write):
    # parse the input data and build locus data structure
    # all inputs should be lists to support multiple inputs
    Preds = []
    n_models = 0
    sources = {"predictions": set(), "evidence": set()}
    for g in natsorted(gene):
        parsed_genes = gff2dict(os.path.abspath(g), genome)
        cleaned = ensure_unique_names(parsed_genes)
        Preds.append(cleaned)
        n_models += len(cleaned)
        for k, v in cleaned.items():
            sources["predictions"].add(v["source"])
    Proteins = {}
    if protein:
        for p in natsorted(protein):
            Proteins = gffevidence2dict(os.path.abspath(p), Proteins)
            for k, v in Proteins.items():
                sources["evidence"].add(v["source"])
    Transcripts = {}
    if transcript:
        for t in natsorted(transcript):
            Transcripts = gffevidence2dict(os.path.abspath(t), Transcripts)
            for k, v in Transcripts.items():
                sources["evidence"].add(v["source"])
    # merge parsed predictions into single dictionary
    log(f"Merging gene predictions from {len(Preds)} source files\n{sources}")
    Genes = {k: v for d in Preds for k, v in d.items()}
    # now build loci from gene models
    loci, n_loci, pseudo = get_loci(Genes)
    pbs = {}
    for x in pseudo:
        if x[1]["source"] not in pbs:
            pbs[x[1]["source"]] = 1
        else:
            pbs[x[1]["source"]] += 1
    log(
        "Parsed {} gene models into {} loci. Dropped {} genes models that were pseudo [labled as such or internal stop codons]\n{}".format(
            n_models, n_loci, len(pseudo), pbs
        )
    )

    # add some evidence
    if len(Proteins) > 0:
        add_evidence(loci, Proteins, source="proteins")
    if len(Transcripts) > 0:
        add_evidence(loci, Transcripts, source="transcripts")
    return loci


def filter4evidence(data, n_genes=3, n_evidence=2):
    # here I want to return loci with at least n_genes and n_evidences
    filt = {}
    n_filt = 0
    for contig, obj in data.items():
        for strand in ["+", "-"]:
            for m in obj[strand]:
                if len(m["genes"]) >= n_genes:
                    if (len(m["transcripts"]) + len(m["proteins"])) >= n_evidence:
                        n_filt += 1
                        if contig not in filt:
                            if strand == "+":
                                filt[contig] = {"+": [m], "-": []}
                            else:
                                filt[contig] = {"-": [m], "+": []}
                        else:
                            if strand not in filt[contig]:
                                filt[contig][strand] = [m]
                            else:
                                filt[contig][strand].append(m)
    return filt, n_filt


def score_by_evidence(locus, weights={}, derived=[]):
    results = {}
    for gene in locus["genes"]:
        score = {"proteins": [], "transcripts": []}
        name, source, coords, cstart = gene
        if source not in derived:
            for s in ["proteins", "transcripts"]:
                for x in locus[s]:
                    q_name, q_source, q_coords = x
                    score[s].append(
                        score_evidence(
                            coords[0], q_coords, weight=weights.get(source, 1)
                        )
                    )
            results[name] = {
                "protein_evidence_score": sum(score["proteins"]),
                "transcript_evidence_score": sum(score["transcripts"]),
            }
        else:
            results[name] = {
                "protein_evidence_score": 0.0,
                "transcript_evidence_score": 0.0,
            }
    return results


def order_sources(locus):
    results = {}
    # print(locus['genes'])
    # first lets check if all predictions are the same
    allcoords = [x[2][0] for x in locus["genes"]]
    if len([list(i) for i in set(tuple(i) for i in allcoords)]) > 0:
        for gene in locus["genes"]:
            score = {"proteins": [], "transcripts": []}
            name, source, coords, cstart = gene
            for s in ["proteins", "transcripts"]:
                for x in locus[s]:
                    q_name, q_source, q_coords = x
                    score[s].append(score_evidence(coords[0], q_coords, weight=1))
            results[name] = {
                "source": source,
                "coords": coords[0],
                "score": sum(score["proteins"]) + sum(score["transcripts"]),
            }
    else:
        results = {
            locus["genes"][0][0]: {
                "source": locus["genes"][0][1],
                "coords": locus["genes"][0][2],
                "score": 0,
            }
        }
    return results


def debug_gff_writer(outfile, contig, strand, data):
    # input is a sorted list of gene models from best_model function
    with zopen(outfile, mode="a") as out:
        for x in data:
            name, d = x
            if strand == "+":
                sortedCoords = sorted(d["coords"], key=lambda tup: tup[0])
            else:
                sortedCoords = sorted(d["coords"], key=lambda tup: tup[0], reverse=True)
            g_start = min(min(sortedCoords))
            g_end = max(max(sortedCoords))
            out.write(
                "{}\t{}\tgene\t{}\t{}\t{:.4f}\t{}\t.\tID={};Note={},{},reasonable:{};\n".format(
                    contig,
                    d["source"],
                    g_start,
                    g_end,
                    d["score"],
                    strand,
                    name,
                    d["locus"],
                    d["all_scores"],
                    d["check"],
                )
            )
            out.write(
                "{}\t{}\tmRNA\t{}\t{}\t.\t{}\t.\tID={}.mrna;Parent={};Note={},{},reasonable:{};\n".format(
                    contig,
                    d["source"],
                    g_start,
                    g_end,
                    strand,
                    name,
                    name,
                    d["locus"],
                    d["all_scores"],
                    d["check"],
                )
            )
            # now write exon and CDS features
            num_exons = len(sortedCoords)
            for x in range(0, num_exons):
                ex_num = x + 1
                out.write(
                    "{:}\t{:}\texon\t{:}\t{:}\t.\t{:}\t.\tID={:}.exon{:};Parent={:}.mrna;\n".format(
                        contig,
                        d["source"],
                        sortedCoords[x][0],
                        sortedCoords[x][1],
                        strand,
                        name,
                        ex_num,
                        name,
                    )
                )
            current_phase = d["codon_start"] - 1
            for y in range(0, num_exons):
                out.write(
                    "{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds;Parent={:}.mrna;\n".format(
                        contig,
                        d["source"],
                        sortedCoords[y][0],
                        sortedCoords[y][1],
                        strand,
                        current_phase,
                        name,
                        name,
                    )
                )
                current_phase = (
                    current_phase
                    - (int(sortedCoords[y][1]) - int(sortedCoords[y][0]) + 1)
                ) % 3
                if current_phase == 3:
                    current_phase = 0


def reasonable_model(
    coords, min_protein=30, min_exon=3, min_intron=10, max_intron=-1, max_exon=-1
):
    # look at coordinates and check some stats
    exon_lengths = []
    intron_lengths = []
    sortedCoords = sorted(coords, key=lambda tup: tup[0])
    for i, x in enumerate(sortedCoords):
        exon_lengths.append(x[1] - x[0])
        if i > 0:
            intron_lengths.append(x[0] - sortedCoords[i - 1][1])
    if min(exon_lengths) < min_exon:
        return "exon < {}".format(min_exon)
    if max_exon > 0 and max(exon_lengths) > max_exon:
        return "exon > {}".format(max_exon)
    if len(intron_lengths) > 0:
        if min(intron_lengths) < min_intron:
            return "intron < {}".format(min_intron)
        if max_intron > 0 and max(intron_lengths) > max_intron:
            return "intron > {}".format(max_intron)
    if (sum(exon_lengths) / 3) < min_protein:
        return "protein is < {} amino acids".format(min_protein)
    return True


def refine_cluster(locus, derived=[]):
    # get sources here, but ignore models derived from evidence, ie we
    # just want to find ab initio models here
    sources = [x[1] for x in locus["genes"] if x[1] not in derived]
    if len(sources) > 0:
        src, n = Counter(sources).most_common(1)[0]
        if n > 1:
            # [['contig_76-snap.4', 'snap', [[(11710, 11713), (12068, 12155), (12543, 12881), (15021, 15918)]], 1]]
            seeds = []
            queries = []
            locs = []
            for x in locus["genes"]:
                if x[1] == src:
                    coords = (min(x[2][0])[0], max(x[2][0])[1])
                    if len(locs) > 0:
                        unique = True
                        for z in locs:
                            if get_overlap(coords, z) > 0:
                                unique = False
                        if unique:
                            locs.append(coords)
                            if isinstance(x, list):
                                seeds.append(x)
                            else:
                                seeds.append([x])
                        else:
                            queries.append(x)
                    else:
                        locs.append(coords)
                        if isinstance(x, list):
                            seeds.append(x)
                        else:
                            seeds.append([x])
                else:
                    queries.append(x)
            if len(seeds) > 1:
                # get final data strucutre
                final = {}
                for x, y in enumerate(seeds):
                    final[x] = [y]
                # if the seeds overlap, then forget it, its a single locus
                # now we can assign the queries to seeds
                for y in queries:
                    for idx, w in enumerate(locs):
                        o = get_overlap([min(y[2][0])[0], max(y[2][0])[1]], w)
                        if o > 0:
                            if isinstance(y, list):
                                final[idx].append(y)
                            else:
                                final[idx].append([y])
                return final
            else:
                return False
        else:
            return False
    else:
        return False


def best_model(
    locus_name,
    contig,
    strand,
    locus,
    debug=False,
    weights={},
    order={},
    min_exon=3,
    min_intron=10,
    max_intron=-1,
    max_exon=-1,
    evidence_derived_models=[],
):
    # get evidence scores
    evidence_scores = score_by_evidence(locus, derived=evidence_derived_models)
    # get aed scores
    de_novo_aed_scores = de_novo_distance(locus)
    # check if could be multiple
    sub_loci = refine_cluster(locus, derived=evidence_derived_models)
    if sub_loci:
        # print("---------------")
        # print(f"Sub clusters identified: {len(sub_loci)}")
        r = []
        seen = set()
        non_redundant = []
        for i, z in sub_loci.items():
            # print(f"cluster {i}: {z}")
            tmp_locus = {
                "genes": z,
                "locus": locus["locus"],
                "proteins": locus["proteins"],
                "transcripts": locus["transcripts"],
                "repeats": locus["repeats"],
            }
            tmp_results = score_aggregator(
                locus_name,
                tmp_locus,
                weights,
                order,
                de_novo_aed_scores,
                evidence_scores,
                min_exon=min_exon,
                min_intron=min_intron,
                max_intron=max_intron,
                max_exon=max_exon,
            )
            tmp_sorted = sorted(
                tmp_results.items(), key=lambda e: e[1]["score"], reverse=True
            )
            for model in tmp_sorted:
                if model[0] not in seen:
                    seen.add(model[0])
                    non_redundant.append(model)
            # drop any models that aren't reasonable, even if highest scoring
            tmp_filtered = [x for x in tmp_sorted if x[1]["check"] is True]
            if len(tmp_filtered) > 0:
                r.append(tmp_filtered)
        # now we need to go through the results and choose wisely, single or multiple models
        if len(r) > 0:
            try:
                besties = set(
                    [
                        (
                            x[0][0],
                            (
                                min(x[0][1]["coords"])[0],
                                max(x[0][1]["coords"])[1],
                            ),
                            max(x[0][1]["coords"])[1] - min(x[0][1]["coords"])[0],
                            i,
                        )
                        for i, x in enumerate(r)
                    ]
                )
            except IndexError:
                print(r)
                raise SystemExit(1)
            # print(r)
            # print("best models:")
            # print(besties)
            # simple first, if best hits all same model
            if len(besties) == 1:  # then just single model to return
                bm = [r[0][0]]
            elif len(besties) > 1:  # more work here
                bms = solve_sub_loci(besties)
                bm = []
                for b in bms:
                    bm.append(r[b[3]][0])
        else:
            bm = []
    else:  # these are singles, shortcut the sub cluster routine
        # now go through genes and combine scores
        results = score_aggregator(
            locus_name,
            locus,
            weights,
            order,
            de_novo_aed_scores,
            evidence_scores,
            min_exon=min_exon,
            min_intron=min_intron,
            max_intron=max_intron,
            max_exon=max_exon,
        )
        best_result = sorted(results.items(), key=lambda e: e[1]["score"], reverse=True)
        non_redundant = best_result
        best_result_filtered = [x for x in best_result if x[1]["check"] is True]
        # check if we need to break any ties
        if len(best_result_filtered) > 0:
            best_score = best_result_filtered[0][1]["score"]
            anyties = [x for x in best_result_filtered if x[1]["score"] >= best_score]
            if len(anyties) > 1:
                breaktie = sorted(
                    anyties, key=lambda x: order[x[1]["source"]], reverse=True
                )
                bm = [random.choice(breaktie)]
            else:
                bm = [anyties[0]]
        else:
            bm = []
    # if debug than write GFF3 to file
    if debug:
        debug_gff_writer(debug, contig, strand, non_redundant)
    return bm


def solve_sub_loci(result):
    # given a set of tuples (id, coords, length, index)
    # {('contig_1-g4', (10518, 12217), 1699, 0), ('contig_1-snap.5', (10801, 16904), 6103, 1)}
    # need to return non-overlapping coords to fill the locus space
    keep = set()
    sbl = sorted(list(result), key=lambda x: x[2], reverse=True)
    keep.add(sbl[0])
    for x in sbl[1:]:
        k_coords = list(list(zip(*keep))[1])
        ovlps = []
        for y in k_coords:
            ovlps.append(get_overlap(x[1], y))
        if len(set(ovlps)) == 1 and ovlps[0] == 0:
            keep.add(x)
    return list(keep)


def score_aggregator(
    locus_name,
    locus,
    weights,
    order,
    de_novo_aed_scores,
    evidence_scores,
    min_exon=3,
    min_intron=10,
    max_intron=-1,
    max_exon=-1,
):
    results = {}
    for gene in locus["genes"]:
        try:
            name, source, coords, cstart = gene
        except ValueError:
            print("ERROR parsing gene object")
            print(gene)
            raise SystemExit(1)
        # check if gene is reasonable
        check = reasonable_model(
            coords[0],
            min_exon=min_exon,
            min_intron=min_intron,
            max_intron=max_intron,
            max_exon=max_exon,
        )
        user_score = weights.get(source, 1) * 6
        order_score = order.get(source, 0)
        aed_score = de_novo_aed_scores.get(name, 0) * 2
        evi = evidence_scores.get(name)
        prot_score = evi["protein_evidence_score"]
        trans_score = evi["transcript_evidence_score"]
        # if user_score is 0, then make total score 0 as user has indicated not to use those models
        if user_score == 0:
            total_score = 0
        else:
            total_score = user_score + aed_score + prot_score + trans_score
        d = {
            "source": source,
            "coords": coords[0],
            "check": check,
            "score": total_score,
            "locus": locus_name,
            "codon_start": cstart,
            "all_scores": "total:{:.2f},user:{},aed:{:.4f},prot:{},trans:{},order:{}".format(
                total_score, user_score, aed_score, prot_score, trans_score, order_score
            ),
        }
        results[name] = d
    return results


def cluster_by_aed(locus, score=0.005):
    distances = calculate_gene_distance(locus)
    # try to cluster, dict is ordered so run through in order
    clusters = {}  # start with identical models
    seen = {}
    for k, v in distances.items():
        # print(k, v)
        if k in seen:
            centroid = seen.get(k)
        else:
            centroid = k
        for g, a in v.items():
            if g in seen:
                centroid = seen.get(g)
            if a <= score:
                if centroid == g:
                    continue
                if centroid not in clusters:
                    clusters[centroid] = [g]
                else:
                    if g not in clusters[centroid]:
                        clusters[centroid].append(g)
                seen[g] = centroid
                if g in clusters:
                    del clusters[g]
            else:
                if g not in clusters and g not in seen:
                    clusters[g] = []
    # print(clusters)
    final = []
    accounted = 0
    for k, v in sorted(clusters.items(), key=lambda x: len(x[1]), reverse=True):
        all = [k] + v
        accounted += len(all)
        final.append(all)
    assert len(locus["genes"]) == accounted
    return final


def map_coords(g_coords, e_coords):
    # map the e_coords onto the g_coords
    r = [
        [],
    ] * len(g_coords)
    refInterlap = interlap.InterLap(g_coords)
    for e in e_coords:
        if e in refInterlap:
            hit = list(refInterlap.find(e))
            if len(hit) > 1:
                continue
            i = g_coords.index(hit[0])
            diff = np.subtract(e, hit[0])
            r[i] = list(diff)
    return r


def score_evidence(g_coords, e_coords, weight=2):
    # simple score: 10 to 0
    # 1 == e_coords contained and match intron/exons
    # 0.5 == e_coords partially contained
    # 0 == e_coords not contained or contained but intron/exon boundaries do not match
    # coords are list of tuples
    # g_coords is single gene model
    # e_coords is single evidence coordinates
    if g_coords == e_coords:
        return 10 * weight
    emap = map_coords(g_coords, e_coords)
    score = 0
    if len(g_coords) > 1:
        mult = []
        for x in emap:
            if not x:
                mult.append(0)
            else:
                if x[0] >= 0 and x[1] <= 0:
                    mult.append(10)
                else:
                    mult.append(5)
        score = (sum(mult) // len(g_coords)) * weight
    else:  # single exon gene
        if emap[0]:
            if emap[0][0] >= 0 and emap[0][1] <= 0:
                score = 10 * weight
            else:
                score = 5 * weight
        else:
            score = 0
    return score


def calculate_source_order(data):
    # filter data for evidence and calculate the global accuracy of the sources
    # this will output a rank order of sources to use when there is no evidence
    f_data, n_filt = filter4evidence(data)
    if n_filt > 0:
        cweights = {}
        for contig, obj in f_data.items():
            for strand in ["+", "-"]:
                for locus in obj[strand]:
                    all = order_sources(locus)
                    for k, v in all.items():
                        if v["source"] not in cweights:
                            cweights[v["source"]] = [v["score"]]
                        else:
                            cweights[v["source"]].append(v["score"])
        sweights = []
        avgweights = []
        for k, v in cweights.items():
            aweight = sum(v) // len(v)
            sweights.append((k, aweight))
            avgweights.append(aweight)
        min_weight = min(avgweights)
        order = (
            OrderedDict()
        )  # y[0] for y in sorted(sweights.items(), key=lambda x: x[1], reverse=True)]
        for z in sorted(sweights, key=lambda x: x[1], reverse=True):
            order[z[0]] = z[1]  # round(z[1] / min_weight)
    else:  # there are no loci with evidence, so return all the sources instead and set all to 1
        order = OrderedDict()
        for contig, obj in data.items():
            for strand in ["+", "-"]:
                for locus in obj[strand]:
                    for gene in locus["genes"]:
                        name, source, coords, cstart = gene
                        if source not in order:
                            order[source] = 1
    return order, n_filt


def calculate_gene_distance(locus):
    # caclulate AED for all pairs of gene predictions
    dist = {}
    for y in locus["genes"]:
        dist[y[0]] = {}
    for x in itertools.combinations(locus["genes"], 2):
        a, b = x
        a_name, a_source, a_coord, _ = a
        b_name, b_source, b_coord, _ = b
        aed = getAED(b_coord[0], a_coord[0])
        dist[a_name][b_name] = aed
        dist[b_name][a_name] = aed
    return dist


def src_scaling_factor(obj):
    # from the aed de novo distance, check how many other sources overlap with
    src = []
    for k, v in obj.items():
        if v < 1.0:  # means some overlap
            src.append(
                k.rsplit(".", 1)[-1]
            )  # this works because src slug is appended to each name
    dups = Counter(src)
    if len(dups) > 0:
        # there are some legit cases where this is perhaps problematic....
        # but we will return here a fraction of the srcs that are 1 and use as scaling factor
        ones = [x for x in dups.values() if x == 1]
        return len(ones) / len(dups)
    else:
        return 1.0


def de_novo_distance(locus):
    results = {}
    if len(locus["genes"]) > 1:
        # calculate a score for each gene based on AED similaritiy with other models
        distances = calculate_gene_distance(locus)
        # add up the aeds, subtract total aed from total gene models
        # idea here is that higher value is a model that has similarity with other predictions
        # this is not as fine-grained as EVM's exon and splice site overlap method, where each exon/intron is evaluated
        for k, v in distances.items():
            total = float(sum(v.values()))
            # we can measure the sources of overlap via the unique slug in each value
            scaling_factor = src_scaling_factor(v)
            assert (
                len(locus["genes"]) == len(v) + 1
            ), "ERROR, AED distance calculation failed"
            results[k] = ((len(locus["genes"])) - total) * scaling_factor
    else:  # single gene, so value here is 0
        geneID = locus["genes"][0][0]
        results[geneID] = 0.0
    return results


def getAED(query, reference):
    """
    function to calcuate annotation edit distance between two transcript coordinates
    AED = 1 - (SN + SP / 2)
    SN = fraction of ref predicted
    SP = fraction prediction overlapping the ref
    """

    def _length(listTup):
        len = 0
        for i in listTup:
            len += abs(i[0] - i[1])
        return len

    # check if identical
    if query == reference:
        return 0.000
    # make sure sorted
    rLen = _length(reference)
    qLen = _length(query)
    refInterlap = interlap.InterLap(reference)
    FPbases = 0
    FNbases = 0
    TPbases = 0
    refSeen = []
    refPerfect = []
    for exon in query:
        if exon in refInterlap:  # exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                refSeen.append(h)
                # will return array of exon minus hit at each pos
                diff = np.subtract(exon, h)
                if diff[0] == 0 and diff[1] == 0:  # then exon is perfect match
                    TPbases += abs(h[0] - h[1])
                    refPerfect.append(h)
                elif diff[0] <= 0 and diff[1] >= 0:  # then query exon covers ref exon
                    cov = abs(h[0] - h[1])
                    FPbases += abs(diff[0])
                    FPbases += abs(diff[1])
                    TPbases += cov
                elif diff[0] <= 0 and diff[1] < 0:  # means query partial covers ref
                    cov = abs(h[0] - exon[1])
                    FPbases += abs(diff[0])
                    FNbases += abs(diff[1])
                    TPbases += cov
                elif diff[0] > 0 and diff[1] >= 0:  # means query partial covers ref
                    cov = abs(exon[0] - h[1])
                    FPbases += abs(diff[1])
                    FNbases += abs(diff[0])
                    TPbases += cov
                elif diff[0] > 0 and diff[1] < 1:
                    cov = abs(exon[0] - exon[1])
                    FNbases += abs(diff[1])
                    FNbases += abs(diff[0])
                    TPbases += cov
        else:
            FPbases += abs(exon[0] - exon[1])
    # last thing is to double check reference for missed exons
    for ex in refInterlap:
        if ex not in refSeen:
            FNbases += abs(ex[0] - ex[1])
    # calc AED
    SN = (rLen - FNbases) / rLen
    SP = (qLen - FPbases) / qLen
    try:
        AED = 1 - ((SN + SP) / 2)
    except ZeroDivisionError:
        AED = 1.00
    return AED


def gff_writer(input, output):
    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # open output for writing and write header
    with zopen(output, mode="w") as outfile:
        outfile.write("##gff-version 3\n")
        # sort the annotations by contig and start location
        sGenes = natsorted(iter(input.items()), key=_sortDict)
        sortedGenes = OrderedDict(sGenes)
        counter = 1
        for k, d in list(sortedGenes.items()):
            locusTag = "CGM_{}".format(str(counter).zfill(6))
            if d["strand"] == "+":
                sortedCoords = sorted(d["coords"], key=lambda tup: tup[0])
            else:
                sortedCoords = sorted(d["coords"], key=lambda tup: tup[0], reverse=True)
            g_start, g_end = d["location"]
            outfile.write(
                "{}\tGFFtk\tgene\t{}\t{}\t{:.4f}\t{}\t.\tID={};Note={} derived from {} model {};\n".format(
                    d["contig"],
                    g_start,
                    g_end,
                    d["score"],
                    d["strand"],
                    locusTag,
                    d["locus"],
                    d["source"],
                    k,
                )
            )
            outfile.write(
                "{}\tGFFtk\tmRNA\t{}\t{}\t.\t{}\t.\tID={}.mrna;Parent={};\n".format(
                    d["contig"], g_start, g_end, d["strand"], locusTag, locusTag
                )
            )
            # now write exon and CDS features
            num_exons = len(sortedCoords)
            current_phase = d["codon_start"] - 1
            for x in range(0, num_exons):
                ex_num = x + 1
                outfile.write(
                    "{}\tGFFtk\texon\t{}\t{}\t.\t{}\t.\tID={}.exon{};Parent={}.mrna;\n".format(
                        d["contig"],
                        sortedCoords[x][0],
                        sortedCoords[x][1],
                        d["strand"],
                        locusTag,
                        ex_num,
                        locusTag,
                    )
                )
                outfile.write(
                    "{}\tGFFtk\tCDS\t{}\t{}\t.\t{}\t{}\tID={}.cds;Parent={}.mrna;\n".format(
                        d["contig"],
                        sortedCoords[x][0],
                        sortedCoords[x][1],
                        d["strand"],
                        current_phase,
                        locusTag,
                        locusTag,
                    )
                )
                current_phase = (
                    current_phase
                    - (int(sortedCoords[x][1]) - int(sortedCoords[x][0]) + 1)
                ) % 3
                if current_phase == 3:
                    current_phase = 0
            counter += 1

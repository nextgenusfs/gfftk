import itertools
import json
import multiprocessing
import os
import random
import sys
import uuid
from collections import Counter, OrderedDict, defaultdict

import numpy as np
from natsort import natsorted

from . import interlap
from .fasta import FASTA, fastaparser
from .gff import gff2dict
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
    check_inputs([args.fasta] + args.genes + args.proteins + args.transcripts + [args.repeats])

    # Get number of processes for parallel execution
    num_processes = args.num_processes if hasattr(args, "num_processes") else None

    # Check if using dynamic programming approach
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
        num_processes=num_processes,
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
    min_intron=11,
    max_intron=-1,
    max_exon=-1,
    evidence_derived_models=[],
    num_processes=None,
    log=sys.stderr.write,
):
    """
    Generate consensus gene models from multiple gene prediction sources and evidence.

    This function is the main entry point for the consensus module. It takes gene predictions
    from multiple sources, along with protein and transcript evidence, and generates consensus
    gene models by selecting the best model at each locus based on evidence and source weights.

    The function performs the following steps:
    1. Parse input GFF3 files and cluster gene models into loci
    2. Calculate source weights based on evidence if tiebreakers="calculated"
    3. Select the best gene model at each locus based on evidence and source weights
    4. Filter out gene models that overlap with repeats (if repeats are provided)
    5. Write the consensus gene models to a GFF3 file

    Parameters:
    -----------
    fasta : str
        Path to the genome FASTA file
    genes : list
        List of paths to gene prediction GFF3 files
    proteins : list
        List of paths to protein alignment GFF3 files
    transcripts : list
        List of paths to transcript alignment GFF3 files
    weights : list
        List of source:weight pairs for weighting gene prediction sources
    out : str
        Path to the output GFF3 file
    debug : bool or str, optional
        Whether to print debug information or path to debug GFF file (default: False)
    minscore : bool or int, optional
        Minimum score threshold for gene models, or False to calculate automatically (default: False)
    repeats : bool or str, optional
        Path to repeats GFF3 or BED file, or False to skip repeat filtering (default: False)
    repeat_overlap : int, optional
        Maximum percentage of gene model that can overlap with repeats (default: 90)
    tiebreakers : str, optional
        Method for calculating source weights, either "calculated" or "user" (default: "calculated")
    min_exon : int, optional
        Minimum exon length in nucleotides (default: 3)
    min_intron : int, optional
        Minimum intron length in nucleotides (default: 11)
    max_intron : int, optional
        Maximum intron length in nucleotides, or -1 for no limit (default: -1)
    max_exon : int, optional
        Maximum exon length in nucleotides, or -1 for no limit (default: -1)
    evidence_derived_models : list, optional
        List of sources that are derived from evidence and should be treated differently (default: [])
    num_processes : int or None, optional
        Number of processes to use for parallel execution, or None for sequential (default: None)
    log : callable, optional
        Function to use for logging (default: sys.stderr.write)

    Returns:
    --------
    dict
        Dictionary of consensus gene models, where keys are gene IDs and values are dictionaries
        containing gene model information (contig, location, strand, source, coords, etc.)
    """
    log("GFFtk consensus will generate the best gene model at each locus")
    if weights:
        WEIGHTS = {}
        for x in weights:
            if ":" in x:
                source, w = x.split(":", 1)
            else:
                source = x
                w = 1
            # Try to convert weight to integer, but use default of 1 if it fails
            try:
                WEIGHTS[source] = int(w)
            except (ValueError, TypeError):
                WEIGHTS[source] = 1
                if debug is True:
                    log(
                        f"Warning: Could not convert weight for {source} to integer, using default of 1"
                    )
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
    if debug is True:
        locus_bed = zopen(out + ".loci.bed", mode="w")
    else:
        locus_bed = open(os.devnull, "w")
    split_stats = {}

    # Prepare locus arguments for parallel processing
    locus_args = []
    locus_info = []

    # Create a single FASTA object to be shared across all loci
    fasta_obj = FASTA(fasta)

    for contig, obj in data.items():
        for strand in ["+", "-"]:
            for locus in obj[strand]:
                name = "locus_{}".format(counter)
                counter += 1

                # Store locus info for later use
                locus_info.append((name, contig, strand, locus))

                # Add FASTA object to locus for validation
                locus_with_fasta = locus.copy()
                locus_with_fasta["fasta"] = fasta_obj

                # Prepare arguments for best_model
                locus_args.append(
                    (
                        name,
                        contig,
                        strand,
                        locus_with_fasta,
                        WEIGHTS,
                        order,
                        debug,
                        min_exon,
                        min_intron,
                        max_intron,
                        max_exon,
                        evidence_derived_models,
                    )
                )

    # Process loci in parallel if num_processes is specified
    if num_processes is not None and num_processes > 1:
        log(f"Processing {len(locus_args)} loci using {num_processes} processes")
        with multiprocessing.Pool(processes=num_processes) as pool:
            results = pool.starmap(best_model_default, locus_args)
    else:
        # Process loci sequentially
        results = [best_model_default(*args) for args in locus_args]

    # Process results
    for i, keepers in enumerate(results):
        name, contig, strand, locus = locus_info[i]

        if not keepers:
            continue

        if len(keepers) not in split_stats:
            split_stats[len(keepers)] = 1
        else:
            split_stats[len(keepers)] += 1

        # write locus bed file
        locus_bed.write(
            f"{contig}\t{locus['locus'][0]}\t{locus['locus'][1]}\t{name};n_genes={len(keepers)}\t0\t{strand}\n"
        )
        if contig not in consensus:
            consensus[contig] = {}
        for keep in keepers:
            if keep[0] not in consensus[contig]:
                # Debug: Print the source field
                if debug is True and "source" in keep[1] and keep[1]["source"].startswith("DP:"):
                    sys.stderr.write(f"Found DP model: {keep[0]}, source: {keep[1]['source']}\n")

                consensus[contig][keep[0]] = {
                    "strand": strand,
                    "locus": name,
                    "codon_start": keep[1]["codon_start"],
                    "coords": keep[1]["coords"],
                    "score": keep[1]["score"],
                    "source": keep[1]["source"],
                }
                # Ensure coordinates are integers
                try:
                    # Extract min and max coordinates
                    coords = keep[1]["coords"]

                    # Handle different coordinate formats
                    all_coords = []
                    for coord_group in coords:
                        if isinstance(coord_group, list):
                            for coord in coord_group:
                                if isinstance(coord, tuple) and len(coord) >= 2:
                                    all_coords.append((int(coord[0]), int(coord[1])))
                        elif isinstance(coord_group, tuple) and len(coord_group) >= 2:
                            all_coords.append((int(coord_group[0]), int(coord_group[1])))

                    # Skip if no valid coordinates
                    if not all_coords:
                        # log(f"No valid coordinates for gene {name}\n")
                        continue

                    # Get min and max coordinates
                    min_coord = min(coord[0] for coord in all_coords)
                    max_coord = max(coord[1] for coord in all_coords)

                    if strand == "+":
                        plus_consensus[contig].add(
                            (
                                min_coord,
                                max_coord,
                                name,
                                keep[1]["score"],
                            )
                        )
                    elif strand == "-":
                        minus_consensus[contig].add(
                            (
                                min_coord,
                                max_coord,
                                name,
                                keep[1]["score"],
                            )
                        )
                except (ValueError, TypeError) as e:
                    # Log the error and skip this gene
                    log(f"Error adding gene {name} to consensus: {e}\n")
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
                # Extract coordinates safely
                coords = safe_extract_coordinates(v["coords"])
                if coords is None:
                    # Skip if coordinates can't be extracted
                    continue

                min_coord, max_coord = coords

                if v["strand"] == "+":
                    hits = list(minus_consensus[c].find((min_coord, max_coord)))
                elif v["strand"] == "-":
                    hits = list(plus_consensus[c].find((min_coord, max_coord)))
                if len(hits) > 0:
                    contain = [contained((min_coord, max_coord), (x[0], x[1])) for x in hits]
                    if any(contain):
                        continue
                total += 1
                # Preserve the original source name, including the DP: prefix for composite models
                c_source = v["source"]
                if c_source not in sources:
                    sources[c_source] = 1
                else:
                    sources[c_source] += 1
                v["contig"] = c
                v["name"] = m
                v["location"] = min_coord, max_coord
                filtered[m] = v

    # now we can filter models in repeat regions
    if repeats:
        final, dropped_n = filter_models_repeats(
            fasta, repeats, filtered, filter_threshold=repeat_overlap, log=log
        )
        log("{} gene models were dropped due to repeat overlap".format(dropped_n))
    else:
        final = filtered

    # Check for composite models
    composite_models = {}
    for k, v in list(consensus.items()):
        for m, d in list(v.items()):
            if "source" in d and d["source"].startswith("DP:"):
                if d["score"] > score_threshold:
                    # Extract coordinates safely
                    coords = safe_extract_coordinates(d["coords"])
                    if coords is None:
                        # Skip if coordinates can't be extracted
                        continue

                    min_coord, max_coord = coords

                    # Check for overlapping models on opposite strand
                    if d["strand"] == "+":
                        hits = list(minus_consensus[k].find((min_coord, max_coord)))
                    elif d["strand"] == "-":
                        hits = list(plus_consensus[k].find((min_coord, max_coord)))

                    if len(hits) > 0:
                        contain = [contained((min_coord, max_coord), (x[0], x[1])) for x in hits]
                        if any(contain):
                            continue

                    # Add the composite model to the final set
                    d["contig"] = k
                    d["name"] = m
                    d["location"] = min_coord, max_coord
                    composite_models[m] = d

                    # Update the sources count
                    if d["source"] not in sources:
                        sources[d["source"]] = 1
                    else:
                        sources[d["source"]] += 1

    # Add the composite models to the final set
    final.update(composite_models)

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
    """
    Filter gene models based on their overlap with repeat regions.

    This function filters out gene models that have a significant overlap with repeat regions.
    It builds an interlap object from the repeat file (GFF3 or BED) and calculates the
    percentage of each gene model that overlaps with repeats. Gene models with overlap
    percentage greater than the filter threshold are removed.

    Parameters:
    -----------
    fasta : str
        Path to the genome FASTA file
    repeats : str
        Path to the repeats GFF3 or BED file
    gene_models : dict
        Dictionary of gene models to filter
    filter_threshold : int, optional
        Maximum percentage of gene model that can overlap with repeats (default: 90)
    log : callable or bool, optional
        Function to use for logging, or False to disable logging (default: False)

    Returns:
    --------
    tuple
        A tuple containing:
        - filtered: Dictionary of filtered gene models
        - dropped: Number of gene models that were filtered out
    """
    dropped = 0
    # return fasta length for stats generation
    seq_length = fasta_length(fasta)
    # build interlap object with repeats
    repeat_length = 0
    repeat_inter = defaultdict(interlap.InterLap)
    if os.path.isfile(repeats):
        if repeats.endswith(".bed"):
            repeat_inter, repeat_length = bed2interlap(repeats, inter=repeat_inter)
        else:
            repeat_inter, repeat_length = gff2interlap(repeats, fasta, inter=repeat_inter)
    pct_repeats = repeat_length / float(seq_length)
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
    """
    Calculate the total length of all sequences in a FASTA file.

    This function reads a FASTA file and sums the lengths of all sequences in the file.
    It can be used to determine the total genome size from a genome FASTA file.

    Parameters:
    -----------
    fasta : str
        Path to the FASTA file

    Returns:
    --------
    int
        Total length of all sequences in the FASTA file
    """
    length = 0
    with zopen(fasta) as infile:
        for title, seq in fastaparser(infile):
            length += len(seq)
    return length


def gff2interlap(infile, fasta, inter=False):
    """
    Parse a GFF3 file and construct a scaffold/gene interlap dictionary.

    This function reads a GFF3 file and creates an interlap object containing the genomic
    features defined in the file. The interlap object allows for efficient overlap queries.

    Parameters:
    -----------
    infile : str
        Path to the GFF3 file
    fasta : str
        Path to the genome FASTA file
    inter : dict or bool, optional
        Existing interlap object to update, or False to create a new one (default: False)

    Returns:
    --------
    tuple
        A tuple containing:
        - inter: Dictionary mapping contig names to interlap objects containing features
        - length: Total length of all features in the GFF3 file
    """
    length = 0
    if not inter:
        inter = defaultdict(interlap.InterLap)
    Genes = gff2dict(infile, fasta)
    for k, v in natsorted(list(Genes.items())):
        inter[v["contig"]].add((v["location"][0], v["location"][1]))
        length += v["location"][0] - v["location"][1]
    return inter, length


def bed2interlap(bedfile, inter=False):
    """
    Parse a BED file and construct a scaffold/feature interlap dictionary.

    This function reads a BED file and creates an interlap object containing the genomic
    features defined in the file. The interlap object allows for efficient overlap queries.

    Parameters:
    -----------
    bedfile : str
        Path to the BED file
    inter : dict or bool, optional
        Existing interlap object to update, or False to create a new one (default: False)

    Returns:
    --------
    tuple
        A tuple containing:
        - inter: Dictionary mapping contig names to interlap objects containing features
        - length: Total length of all features in the BED file
    """
    # load interlap object from a bed file
    length = 0
    if not inter:
        inter = defaultdict(interlap.InterLap)
    with zopen(bedfile) as infile:
        for line in infile:
            line = line.strip()
            contig, start, end = line.split("\t")[:3]
            inter[contig].add((int(start), int(end)))
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
    """
    Cluster genomic features using the interlap.reduce function.

    This function takes an interlap object containing genomic features and clusters them
    based on their coordinates. Features that overlap or are adjacent to each other are
    grouped into the same cluster. The function then assigns the original feature data
    to each cluster.

    Parameters:
    -----------
    obj : interlap.InterLap
        Interlap object containing genomic features

    Returns:
    --------
    list
        List of dictionaries, where each dictionary represents a cluster of features
        Each dictionary contains:
        - locus: Tuple of (start, end) coordinates for the cluster
        - genes: List of gene features in the cluster
        - proteins: Empty list for protein evidence (filled later)
        - transcripts: Empty list for transcript evidence (filled later)
        - repeats: Empty list for repeat features (filled later)
    """
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


def sub_cluster(obj):
    """
    Split a cluster of gene models into sub-clusters based on source.

    This function analyzes a cluster of gene models to determine if it contains multiple
    models from the same source. If it does, it splits the cluster into sub-clusters,
    where each sub-cluster contains models that are more likely to belong together based
    on their overlap.

    Parameters:
    -----------
    obj : list
        List of gene model tuples, where each tuple contains:
        (name, source, coords, codon_start)

    Returns:
    --------
    list
        List of lists, where each inner list contains gene models that belong to the
        same sub-cluster
    """
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
    """
    Check if coordinates in a are completely contained within coordinates in b.

    This function determines if the interval defined by coordinates a is completely
    contained within the interval defined by coordinates b. It handles various edge
    cases and ensures that the coordinates are properly formatted as tuples or lists of integers.

    Parameters:
    -----------
    a : tuple or list
        Tuple or list of (start, end) coordinates to check if contained
    b : tuple or list
        Tuple or list of (start, end) coordinates that might contain a

    Returns:
    --------
    bool
        True if a is completely contained within b, False otherwise
    """
    # check if coords in a are completely contained in b
    try:
        # Ensure a and b are converted to integers
        if isinstance(a, (tuple, list)) and len(a) >= 2:
            a_start = int(a[0])
            a_end = int(a[1])
        else:
            return False

        if isinstance(b, (tuple, list)) and len(b) >= 2:
            b_start = int(b[0])
            b_end = int(b[1])
        else:
            return False

        # Check if a is contained in b
        # a is contained in b if a_start >= b_start and a_end <= b_end
        # but not if they are identical
        if a_start >= b_start and a_end <= b_end and not (a_start == b_start and a_end == b_end):
            return True
    except (ValueError, TypeError):
        # If conversion fails, return False
        return False

    return False


def get_overlap(a, b):
    """
    Calculate the overlap between two genomic intervals.

    This function calculates the number of base pairs that overlap between two genomic
    intervals. If the intervals do not overlap, it returns 0.

    Parameters:
    -----------
    a : tuple or list
        Tuple or list of (start, end) coordinates for the first interval
    b : tuple or list
        Tuple or list of (start, end) coordinates for the second interval

    Returns:
    --------
    int
        Number of base pairs that overlap between the two intervals, or 0 if they don't overlap
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_loci(annot_dict):
    """
    Organize gene models into loci based on genomic coordinates and strand.

    This function takes a dictionary of gene models and organizes them into loci based on
    their genomic coordinates and strand. It creates interlap objects for efficient overlap
    queries and clusters overlapping gene models into loci. It also filters out pseudogenes
    and gene models with multiple stop codons.

    Parameters:
    -----------
    annot_dict : dict
        Dictionary of gene models, where keys are gene IDs and values are dictionaries
        containing gene model information (contig, location, strand, source, CDS, etc.)

    Returns:
    --------
    tuple
        A tuple containing:
        - loci: Hierarchical dictionary of loci organized by contig and strand
          {contig: {"+": [locus1, locus2, ...], "-": [locus1, locus2, ...]}}
        - n_loci: Total number of loci
        - pseudo: List of pseudogenes that were filtered out
    """
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
    Parse evidence alignments from a GFF3 file into a dictionary.

    This function reads a GFF3 file containing evidence alignments (proteins or transcripts)
    and converts it into a dictionary mapping alignment IDs to their information. It handles
    multi-exon alignments by combining exons with the same ID into a single entry.

    Parameters:
    -----------
    file : str
        Path to the GFF3 file containing evidence alignments
    Evi : dict
        Existing dictionary to update with new evidence alignments

    Returns:
    --------
    dict
        Dictionary mapping alignment IDs to their information (target, type, source,
        strand, phase, contig, coords, location, score)
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
                except (IndexError, ValueError):
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
    """
    Add evidence alignments to loci based on genomic coordinates.

    This function associates evidence alignments (proteins or transcripts) with gene loci
    based on their genomic coordinates. It builds interlap objects for efficient overlap
    queries and adds the evidence to each locus that it overlaps with.

    Parameters:
    -----------
    loci : dict
        Hierarchical dictionary of loci organized by contig and strand
        {contig: {"+": [locus1, locus2, ...], "-": [locus1, locus2, ...]}}
    evidence : dict
        Dictionary of evidence alignments
    source : str, optional
        Type of evidence being added, either "proteins" or "transcripts" (default: "proteins")

    Returns:
    --------
    None
        The function modifies the loci dictionary in place
    """
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
    """
    Ensure gene model names are unique by appending a UUID slug.

    This function takes a dictionary of gene models and appends a unique identifier
    to each gene model name to ensure there are no name collisions when combining
    gene models from multiple sources.

    Parameters:
    -----------
    genes : dict
        Dictionary of gene models where keys are gene IDs

    Returns:
    --------
    dict
        Dictionary of gene models with unique IDs
    """
    # takes annotation dictionary and appends unique slug
    slug = str(uuid.uuid4())[:8]
    cleaned = {}
    for k, v in genes.items():
        cleaned[f"{k}.{slug}"] = v
    return cleaned


def parse_data(genome, gene, protein, transcript, log=sys.stderr.write):
    """
    Parse input data files and build a locus data structure.

    This function reads gene prediction files, protein alignment files, and transcript
    alignment files, and organizes them into a hierarchical locus data structure. It
    assigns unique identifiers to gene models to avoid name collisions and tracks the
    sources of all predictions and evidence.

    Parameters:
    -----------
    genome : str
        Path to the genome FASTA file
    gene : list
        List of paths to gene prediction GFF3 files
    protein : list
        List of paths to protein alignment GFF3 files
    transcript : list
        List of paths to transcript alignment GFF3 files
    log : callable, optional
        Function to use for logging (default: sys.stderr.write)

    Returns:
    --------
    dict
        Hierarchical dictionary of loci organized by contig and strand
        {contig: {"+": [locus1, locus2, ...], "-": [locus1, locus2, ...]}}
    """
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
    """
    Filter loci to include only those with sufficient gene models and evidence alignments.

    This function filters the loci data structure to include only loci that have at least
    a specified number of gene models and evidence alignments (proteins + transcripts).
    It is used to identify high-confidence loci for calculating source weights.

    Parameters:
    -----------
    data : dict
        Hierarchical dictionary of loci organized by contig and strand
        {contig: {"+": [locus1, locus2, ...], "-": [locus1, locus2, ...]}}
    n_genes : int, optional
        Minimum number of gene models required in a locus (default: 3)
    n_evidence : int, optional
        Minimum number of evidence alignments (proteins + transcripts) required in a locus (default: 2)

    Returns:
    --------
    tuple
        A tuple containing:
        - filt: Filtered dictionary of loci with sufficient gene models and evidence
        - n_filt: Number of loci that passed the filter
    """
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
    """
    Calculate evidence-based scores for gene models in a locus.

    This function evaluates each gene model in a locus by calculating how well it aligns
    with protein and transcript evidence. It assigns scores based on the alignment quality
    and the weights assigned to different gene prediction sources.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models and evidence for a single locus
        Required keys: 'genes', 'proteins', 'transcripts'
    weights : dict, optional
        Dictionary mapping gene model sources to weight values (default: {})
    derived : list, optional
        List of sources that are derived from evidence and should not be scored (default: [])

    Returns:
    --------
    dict
        Dictionary mapping gene model names to evidence scores:
        - protein_evidence_score: Sum of protein evidence scores
        - transcript_evidence_score: Sum of transcript evidence scores
    """
    results = {}
    for gene in locus["genes"]:
        score = {"proteins": [], "transcripts": []}
        name, source, coords, cstart = gene
        if source not in derived:
            for s in ["proteins", "transcripts"]:
                for x in locus[s]:
                    q_name, q_source, q_coords = x
                    score[s].append(
                        score_evidence(coords[0], q_coords, weight=weights.get(source, 1)),
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
    """
    Calculate evidence-based scores for gene models in a locus for source ranking.

    This function evaluates each gene model in a locus by calculating how well it aligns
    with protein and transcript evidence. Unlike score_by_evidence, this function is used
    specifically for ranking gene prediction sources based on their agreement with evidence.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models and evidence for a single locus
        Required keys: 'genes', 'proteins', 'transcripts'

    Returns:
    --------
    dict
        Dictionary mapping gene model names to information:
        - source: Source of the gene model
        - coords: Coordinates of the gene model
        - score: Combined evidence score for the gene model
    """
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


def reasonable_model(coords, min_protein=30, min_exon=3, min_intron=10, max_intron=-1, max_exon=-1):
    """
    Check if a gene model has reasonable exon and intron lengths.

    This function evaluates a gene model to determine if it has reasonable exon and intron
    lengths based on specified thresholds. It checks minimum exon length, maximum exon length,
    minimum intron length, maximum intron length, and minimum protein length.

    Parameters:
    -----------
    coords : list
        List of (start, end) coordinate tuples for the gene model's exons
    min_protein : int, optional
        Minimum protein length in amino acids (default: 30)
    min_exon : int, optional
        Minimum exon length in nucleotides (default: 3)
    min_intron : int, optional
        Minimum intron length in nucleotides (default: 10)
    max_intron : int, optional
        Maximum intron length in nucleotides, or -1 for no limit (default: -1)
    max_exon : int, optional
        Maximum exon length in nucleotides, or -1 for no limit (default: -1)

    Returns:
    --------
    bool or str
        True if the gene model is reasonable, or a string describing the reason it's not reasonable
    """
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
    """
    Identify potential sub-loci within a locus based on non-overlapping gene models.

    This function analyzes a locus to determine if it contains multiple non-overlapping
    gene models from the same source, which might indicate that the locus should be split
    into multiple sub-loci. It focuses on ab initio gene predictors, which typically don't
    predict overlapping models, so multiple models from the same source in a locus suggest
    the presence of multiple genes.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models and evidence for a single locus
        Required keys: 'genes'
    derived : list, optional
        List of sources that are derived from evidence and should be ignored (default: [])

    Returns:
    --------
    dict or bool
        Dictionary mapping sub-locus indices to lists of gene models if sub-loci are found,
        or False if no sub-loci are identified
    """
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


def best_model_default(
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
    """
    Select the best gene model(s) for a locus using evidence-based scoring.

    Parameters:
    -----------
    locus_name : str
        Name of the locus
    contig : str
        Name of the contig
    strand : str
        Strand of the locus ("+" or "-")
    locus : dict
        Dictionary containing gene models and evidence
    debug : bool or str, optional
        Whether to print debug information or path to debug GFF file
    weights : dict, optional
        Dictionary of weights for different sources
    order : dict, optional
        Dictionary of order values for different sources
    min_exon : int, optional
        Minimum exon length
    min_intron : int, optional
        Minimum intron length
    max_intron : int, optional
        Maximum intron length
    max_exon : int, optional
        Maximum exon length
    evidence_derived_models : list, optional
        List of evidence-derived models
    use_dp : bool, optional
        Whether to use dynamic programming approach (not used in this function)
    allow_multiple : bool, optional
        Whether to allow multiple gene models
    min_gap_size : int, optional
        Minimum gap size for splitting loci

    Returns:
    --------
    list
        List of tuples containing gene model IDs and their information
    """
    # Debug information
    if debug is True:
        sys.stderr.write(f"\n=== Processing locus {locus_name} ===\n")
        sys.stderr.write(f"Contig: {contig}, Strand: {strand}\n")
        sys.stderr.write(f"Number of genes: {len(locus.get('genes', []))}\n")
        sys.stderr.write(f"Number of proteins: {len(locus.get('proteins', []))}\n")
        sys.stderr.write(f"Number of transcripts: {len(locus.get('transcripts', []))}\n")

    # Get evidence scores
    evidence_scores = score_by_evidence(locus, derived=evidence_derived_models)
    # Get AED scores (annotation edit distance)
    de_novo_aed_scores = de_novo_distance(locus)

    # check if could be multiple
    sub_loci = refine_cluster(locus, derived=evidence_derived_models)
    if sub_loci:
        r = []
        seen = set()
        non_redundant = []
        for i, z in sub_loci.items():
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
            tmp_sorted = sorted(tmp_results.items(), key=lambda e: e[1]["score"], reverse=True)
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
            except IndexError as e:
                if debug is True:
                    sys.stderr.write(f"Error processing sub-loci: {e}\n")
                    sys.stderr.write(f"{r}\n")
                # Return an empty list instead of exiting the program
                return []

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
                bm = [random.choice(anyties)]
            else:
                bm = [anyties[0]]
        else:
            bm = []
    if debug is True:
        sys.stderr.write(f"Number of consensus models: {len(bm)}\n")
        sys.stderr.write(f"Consensus models:\n{bm}\n")
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
        order_score = order.get(source, 0) if isinstance(order, dict) else 0
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
    """
    Cluster gene models based on their Annotation Edit Distance (AED).

    This function groups gene models that have very similar exon structures (low AED)
    into clusters. It is used to identify gene models that are essentially the same
    prediction but come from different sources, allowing the consensus module to
    select the best representative from each cluster.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models for a single locus
        Required keys: 'genes'
    score : float, optional
        Maximum AED threshold for considering two gene models as part of the same
        cluster (default: 0.005)

    Returns:
    --------
    list
        List of lists, where each inner list contains gene model IDs that belong to
        the same cluster
    """
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
    """
    Map evidence coordinates onto gene model coordinates.

    This function takes evidence coordinates (protein or transcript alignments) and maps them
    onto gene model coordinates. It calculates the offset between each evidence coordinate
    and the corresponding gene model coordinate, which is used to determine how well the
    evidence aligns with the gene model.

    Parameters:
    -----------
    g_coords : list
        List of (start, end) coordinate tuples for the gene model's exons
    e_coords : list
        List of (start, end) coordinate tuples for the evidence alignments

    Returns:
    --------
    list
        List of lists, where each inner list contains the offset between an evidence
        coordinate and the corresponding gene model coordinate. The list has the same
        length as g_coords, with empty lists for gene model coordinates that don't
        have a corresponding evidence coordinate.
    """
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
    """
    Calculate a score for how well evidence aligns with a gene model.

    This function evaluates how well evidence coordinates (protein or transcript alignments)
    match a gene model's exon structure. It considers both the coverage (percentage of the
    gene model covered by evidence) and the matching of intron junctions (splice sites).

    The scoring system ranges from 0 to 10 (before applying the weight multiplier):
    - 10: Perfect match (evidence exactly matches the gene model)
    - 5-9: Partial match (evidence partially covers the gene model or has some matching junctions)
    - 0: No match (evidence does not overlap with the gene model)

    For multi-exon genes, the score is adjusted based on:
    - Base score from exon coverage (0-10 for each exon)
    - Percent coverage of the entire gene model
    - Ratio of matching intron junctions

    Parameters:
    -----------
    g_coords : list
        List of (start, end) coordinate tuples for the gene model's exons
    e_coords : list
        List of (start, end) coordinate tuples for the evidence alignments
    weight : int, optional
        Weight multiplier to apply to the final score (default: 2)

    Returns:
    --------
    int
        Score indicating how well the evidence supports the gene model,
        ranging from 0 (no support) to higher values (strong support)
    """
    # Enhanced scoring: 10 to 0 with percent coverage and intron junction consideration
    # 1 == e_coords contained and match intron/exons with 100% coverage and all junctions match
    # 0.5 == e_coords partially contained or partial coverage
    # 0 == e_coords not contained or contained but intron/exon boundaries do not match
    # coords are list of tuples
    # g_coords is single gene model
    # e_coords is single evidence coordinates

    # Calculate total gene model length
    total_gene_length = sum(end - start for start, end in g_coords)

    # If exact match, return maximum score
    if g_coords == e_coords:
        return 10 * weight

    # Map evidence coordinates to gene coordinates
    emap = map_coords(g_coords, e_coords)
    score = 0

    # Calculate overlap and coverage
    total_overlap = 0

    # Count matching intron junctions
    matching_junctions = 0
    total_junctions = max(0, len(g_coords) - 1)  # Number of introns in gene model

    if len(g_coords) > 1:  # Multi-exon gene
        # Check for matching intron junctions
        if len(e_coords) > 1:  # Only check junctions if evidence has multiple exons
            for i in range(len(g_coords) - 1):
                # End of current exon and start of next exon define an intron junction
                g_junction_end = g_coords[i][1]
                g_junction_start = g_coords[i + 1][0]

                # Check if this junction is supported by evidence
                for j in range(len(e_coords) - 1):
                    e_junction_end = e_coords[j][1]
                    e_junction_start = e_coords[j + 1][0]

                    # If junction positions match exactly, count it as a match
                    if g_junction_end == e_junction_end and g_junction_start == e_junction_start:
                        matching_junctions += 1
                        break

        # Calculate junction match ratio
        junction_match_ratio = matching_junctions / total_junctions if total_junctions > 0 else 0

        # Continue with existing exon overlap calculation
        mult = []
        for i, x in enumerate(emap):
            if not x:  # No overlap for this exon
                mult.append(0)
                continue

            # Calculate overlap for this exon
            exon_length = g_coords[i][1] - g_coords[i][0]

            if x[0] >= 0 and x[1] <= 0:  # Evidence contained within exon
                # Calculate exact overlap
                overlap = exon_length - abs(x[0]) - abs(x[1])
                total_overlap += overlap
                # Perfect match for this exon
                mult.append(10)
            else:  # Partial overlap
                # Calculate partial overlap
                if x[0] < 0:  # Evidence extends before exon
                    start_overlap = g_coords[i][0]
                else:
                    start_overlap = g_coords[i][0] + x[0]

                if x[1] > 0:  # Evidence extends after exon
                    end_overlap = g_coords[i][1]
                else:
                    end_overlap = g_coords[i][1] + x[1]

                overlap = max(0, end_overlap - start_overlap)
                total_overlap += overlap
                mult.append(5)

        # Calculate base score
        base_score = sum(mult) // len(g_coords)

        # Calculate percent coverage
        percent_coverage = total_overlap / total_gene_length if total_gene_length > 0 else 0

        # Adjust score based on coverage AND junction matches
        # Give equal weight to coverage and junction matches
        score = (
            int(base_score * (0.3 + 0.3 * percent_coverage + 0.4 * junction_match_ratio)) * weight
        )

    else:  # Single exon gene
        exon_length = g_coords[0][1] - g_coords[0][0]

        if emap[0]:  # Some overlap
            if emap[0][0] >= 0 and emap[0][1] <= 0:  # Evidence contained within exon
                # Calculate exact overlap
                overlap = exon_length - abs(emap[0][0]) - abs(emap[0][1])
                percent_coverage = overlap / exon_length if exon_length > 0 else 0
                score = int(10 * (0.5 + 0.5 * percent_coverage)) * weight
            else:  # Partial overlap
                # Calculate partial overlap
                if emap[0][0] < 0:  # Evidence extends before exon
                    start_overlap = g_coords[0][0]
                else:
                    start_overlap = g_coords[0][0] + emap[0][0]

                if emap[0][1] > 0:  # Evidence extends after exon
                    end_overlap = g_coords[0][1]
                else:
                    end_overlap = g_coords[0][1] + emap[0][1]

                overlap = max(0, end_overlap - start_overlap)
                percent_coverage = overlap / exon_length if exon_length > 0 else 0
                score = int(5 * (0.5 + 0.5 * percent_coverage)) * weight
        else:  # No overlap
            score = 0
    return score


def calculate_source_order(data):
    """
    Calculate a rank order of gene prediction sources based on evidence agreement.

    This function analyzes how well each gene prediction source agrees with protein and
    transcript evidence across all loci. It filters the data to include only loci with
    sufficient evidence, calculates scores for each source based on evidence agreement,
    and returns a rank-ordered dictionary of sources with their scores.

    The rank order is used to prioritize gene models when evidence is not available or
    is inconclusive. Sources that generally have better agreement with evidence receive
    higher scores and are ranked higher.

    Parameters:
    -----------
    data : dict
        Hierarchical dictionary of loci organized by contig and strand
        {contig: {"+": [locus1, locus2, ...], "-": [locus1, locus2, ...]}}

    Returns:
    --------
    tuple
        A tuple containing:
        - order: OrderedDict mapping source names to their scores, ordered by score (highest first)
        - n_filt: Number of loci that passed the evidence filter
    """
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
        # min_weight = min(avgweights)  # Not used
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
    """
    Calculate Annotation Edit Distance (AED) between all pairs of gene models in a locus.

    This function computes the AED between each pair of gene models in a locus, which
    measures how similar their exon structures are. The AED scores are used to determine
    which gene models have the most agreement with other models.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models for a single locus
        Required keys: 'genes'

    Returns:
    --------
    dict
        Nested dictionary mapping gene model names to dictionaries of AED scores with other models
        {gene1: {gene2: aed_score, gene3: aed_score, ...}, gene2: {...}, ...}
    """
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
    """
    Calculate a scaling factor based on the diversity of gene prediction sources that agree.

    This function analyzes the AED scores between gene models to determine how many different
    gene prediction sources agree with each other. It returns a scaling factor that reflects
    the proportion of unique sources that have only one model in agreement with others.

    The scaling factor is used to adjust de novo distance scores to favor gene models that
    have agreement across multiple different sources rather than multiple models from the
    same source.

    Parameters:
    -----------
    obj : dict
        Dictionary mapping gene model names to their AED scores with other gene models

    Returns:
    --------
    float
        Scaling factor between 0 and 1, where 1 indicates all sources have only one model
        in agreement with others, and lower values indicate multiple models from the same
        source agree with others
    """
    # from the aed de novo distance, check how many other sources overlap with
    src = []
    for k, v in obj.items():
        if v < 1.0:  # means some overlap
            src.append(k.rsplit(".", 1)[-1])  # this works because src slug is appended to each name
    dups = Counter(src)
    if len(dups) > 0:
        # there are some legit cases where this is perhaps problematic....
        # but we will return here a fraction of the srcs that are 1 and use as scaling factor
        ones = [x for x in dups.values() if x == 1]
        return len(ones) / len(dups)
    else:
        return 1.0


def de_novo_distance(locus):
    """
    Calculate a score for each gene model based on its similarity to other models.

    This function evaluates each gene model in a locus by calculating its Annotation Edit
    Distance (AED) with all other gene models in the locus. It then computes a score that
    reflects how similar the gene model is to other models, with higher scores indicating
    greater similarity.

    This score is used as a proxy for prediction confidence when protein or transcript
    evidence is not available. Gene models that have more agreement with other models
    receive higher scores.

    Parameters:
    -----------
    locus : dict
        Dictionary containing gene models for a single locus
        Required keys: 'genes'

    Returns:
    --------
    dict
        Dictionary mapping gene model names to their de novo distance scores
    """
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
            assert len(locus["genes"]) == len(v) + 1, "ERROR, AED distance calculation failed"
            results[k] = ((len(locus["genes"])) - total) * scaling_factor
    else:  # single gene, so value here is 0
        geneID = locus["genes"][0][0]
        results[geneID] = 0.0
    return results


def getAED(query, reference):
    """
    Calculate Annotation Edit Distance (AED) between two transcript coordinates.

    AED measures the similarity between two gene models by comparing their exon structures.
    It is calculated as 1 - (SN + SP) / 2, where:
    - SN (Sensitivity) is the fraction of reference bases that are correctly predicted
    - SP (Specificity) is the fraction of prediction bases that overlap with the reference

    An AED of 0 means the gene models are identical, while an AED of 1 means they are completely different.

    Parameters:
    -----------
    query : list
        List of (start, end) coordinate tuples for the query gene model's exons
    reference : list
        List of (start, end) coordinate tuples for the reference gene model's exons

    Returns:
    --------
    float
        AED score between 0 (identical) and 1 (completely different)
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
    """
    Write consensus gene models to a GFF3 file.

    This function takes a dictionary of consensus gene models and writes them to a GFF3 file.
    It sorts the gene models by contig and start location, and assigns sequential locus tags
    to each gene model. It also handles the conversion of gene model coordinates to GFF3
    features (gene, mRNA, exon, CDS).

    Parameters:
    -----------
    input : dict
        Dictionary of consensus gene models, where keys are gene IDs and values are dictionaries
        containing gene model information (contig, location, strand, source, coords, etc.)
    output : str
        Path to the output GFF3 file

    Returns:
    --------
    None
        The function writes to the specified output file but does not return a value
    """

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
            # Check if location is valid
            try:
                g_start, g_end = d["location"]
            except (ValueError, TypeError):
                # Skip this gene if location is invalid
                continue

            # Preserve the original source name, including the DP: prefix for composite models
            source = d["source"]
            outfile.write(
                "{}\tGFFtk\tgene\t{}\t{}\t{:.4f}\t{}\t.\tID={};Note={} derived from {} model {};\n".format(
                    d["contig"],
                    g_start,
                    g_end,
                    d["score"],
                    d["strand"],
                    locusTag,
                    d["locus"],
                    source,
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

            # Check if sortedCoords is valid
            valid_coords = []
            for coord in sortedCoords:
                if isinstance(coord, (list, tuple)) and len(coord) >= 2:
                    try:
                        start = int(coord[0])
                        end = int(coord[1])
                        valid_coords.append((start, end))
                    except (ValueError, TypeError):
                        pass

            # Use valid_coords instead of sortedCoords
            sortedCoords = valid_coords
            num_exons = len(sortedCoords)

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
                    current_phase - (int(sortedCoords[x][1]) - int(sortedCoords[x][0]) + 1)
                ) % 3
                if current_phase == 3:
                    current_phase = 0
            counter += 1


def safe_extract_coordinates(coords):
    """
    Safely extract min and max coordinates from a nested coordinate structure.

    This function handles various coordinate formats and structures, extracting the
    minimum and maximum coordinates while gracefully handling errors. It's designed
    to work with potentially malformed or inconsistent coordinate data.

    Parameters:
    -----------
    coords : list or tuple
        Nested coordinate structure (list of lists, tuples, etc.)

    Returns:
    --------
    tuple or None
        A tuple containing (min_coord, max_coord) if extraction succeeds,
        or None if extraction fails
    """
    try:
        # Handle different coordinate formats
        all_coords = []
        for coord_group in coords:
            if isinstance(coord_group, list):
                for coord in coord_group:
                    if isinstance(coord, tuple) or isinstance(coord, list):
                        if len(coord) >= 2:
                            try:
                                all_coords.append((int(coord[0]), int(coord[1])))
                            except (ValueError, TypeError):
                                pass
            elif isinstance(coord_group, tuple) or isinstance(coord_group, list):
                if len(coord_group) >= 2:
                    try:
                        all_coords.append((int(coord_group[0]), int(coord_group[1])))
                    except (ValueError, TypeError):
                        pass

        # Skip if no valid coordinates
        if not all_coords:
            return None

        # Get min and max coordinates
        min_coord = min(coord[0] for coord in all_coords)
        max_coord = max(coord[1] for coord in all_coords)

        return (min_coord, max_coord)
    except Exception:
        return None

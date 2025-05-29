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
    utrs=True,  # New parameter to control UTR extension
    min_utr_length=10,
    max_utr_length=2000,
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
                        # Handle different gene tuple formats
                        if len(gene) == 4:
                            name, source, coords, cstart = gene
                        elif len(gene) == 5:
                            name, source, coords, cstart, full_location = gene
                        elif len(gene) == 7:
                            (
                                name,
                                source,
                                coords,
                                cstart,
                                full_location,
                                utr_5,
                                utr_3,
                            ) = gene
                        elif len(gene) == 9:
                            (
                                start,
                                end,
                                name,
                                source,
                                coords,
                                cstart,
                                full_location,
                                utr_5,
                                utr_3,
                            ) = gene
                        else:
                            raise ValueError(f"Unexpected gene tuple length: {len(gene)}")

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

    # Process results, build locus map file
    locusMap = OrderedDict()
    for i, keepers in enumerate(results):
        name, contig, strand, locus = locus_info[i]
        locusMap[name] = {
            "contig": contig,
            "strand": strand,
            "input_models": locus,
            "consensus_models": [],
        }

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
                consensus[contig][keep[0]] = {
                    "strand": strand,
                    "locus": name,
                    "codon_start": keep[1]["codon_start"],
                    "coords": keep[1]["coords"],
                    "score": keep[1]["score"],
                    "source": keep[1]["source"],
                }
                # Add full_location if available
                if "full_location" in keep[1]:
                    consensus[contig][keep[0]]["full_location"] = keep[1]["full_location"]
                # Add UTR information if available
                if "5UTR" in keep[1]:
                    consensus[contig][keep[0]]["5UTR"] = keep[1]["5UTR"]
                if "3UTR" in keep[1]:
                    consensus[contig][keep[0]]["3UTR"] = keep[1]["3UTR"]
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
                # Extract coordinates safely from CDS for overlap checking
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
                # Preserve the original source name
                c_source = v["source"]
                if c_source not in sources:
                    sources[c_source] = 1
                else:
                    sources[c_source] += 1
                v["contig"] = c
                v["name"] = m
                # Use full gene location (including UTRs) if available, otherwise use CDS coordinates
                if "full_location" in v and v["full_location"] is not None:
                    v["location"] = v["full_location"]
                else:
                    # Fallback to CDS coordinates if no UTR information is available
                    v["location"] = min_coord, max_coord
                # Always store CDS coordinates separately for reference
                v["cds_location"] = min_coord, max_coord

                # Transfer UTR information if it exists
                if "5UTR" in v and v["5UTR"]:
                    # Convert 5UTR list of lists to flat list of tuples
                    five_prime_utrs = []
                    for transcript_utrs in v["5UTR"]:
                        if transcript_utrs:  # Skip empty lists
                            five_prime_utrs.extend(transcript_utrs)
                    if five_prime_utrs:
                        v["five_prime_utr"] = five_prime_utrs

                if "3UTR" in v and v["3UTR"]:
                    # Convert 3UTR list of lists to flat list of tuples
                    three_prime_utrs = []
                    for transcript_utrs in v["3UTR"]:
                        if transcript_utrs:  # Skip empty lists
                            three_prime_utrs.extend(transcript_utrs)
                    if three_prime_utrs:
                        v["three_prime_utr"] = three_prime_utrs

                filtered[m] = v
                # Add the model to the locus map
                locusMap[v["locus"]]["consensus_models"].append(v)

    # output the locusMap
    locusMapOut = out + ".loci.json"
    with open(locusMapOut, "w") as locusMapOut:
        json.dump(locusMap, locusMapOut, indent=2)

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
    # check if utr mode is enabled
    if utrs and transcripts:
        log("Extending gene models with UTRs based on transcript evidence")
        # Parse transcript alignments if not already done
        transcript_alignments = {}
        for t in natsorted(transcripts):
            transcript_alignments = gffevidence2dict(os.path.abspath(t), transcript_alignments)

        # Extend UTRs based on transcript evidence
        final = extend_utrs(
            final,
            transcript_alignments,
            fasta,
            min_utr_length=min_utr_length,
            max_utr_length=max_utr_length,
            log=log,
        )
    else:
        log("UTR extension not enabled, using CDS-only gene models")
        # Copy coords to cds for consistency
        for model_id, model in final.items():
            if "coords" in model:
                model["cds"] = model["coords"].copy()

    # finally we can write to GFF3 format
    gff_writer(final, out)
    log("GFFtk consensus is finished: {}".format(out))
    return final


def check_intron_compatibility(model_coords, transcript_coords, strand):
    """
    Check if transcript has compatible intron/exon boundaries with the gene model.

    Parameters:
    -----------
    model_coords : list
        List of (start, end) tuples for the gene model exons
    transcript_coords : list
        List of (start, end) tuples for the transcript alignment
    strand : str
        Strand of the gene model ('+' or '-')

    Returns:
    --------
    bool
        True if the transcript has compatible intron/exon boundaries, False otherwise
    """
    # Extract intron coordinates from model
    model_introns = []
    sorted_model = sorted(model_coords, key=lambda x: x[0])
    for i in range(len(sorted_model) - 1):
        model_introns.append((sorted_model[i][1], sorted_model[i + 1][0]))

    # If model has no introns, any transcript is compatible
    if not model_introns:
        return True

    # Extract intron coordinates from transcript
    transcript_introns = []
    sorted_transcript = sorted(transcript_coords, key=lambda x: x[0])
    for i in range(len(sorted_transcript) - 1):
        transcript_introns.append((sorted_transcript[i][1], sorted_transcript[i + 1][0]))

    # If transcript has no introns but model does, it's not compatible
    if not transcript_introns:
        return False

    # Check if transcript introns are compatible with model introns
    # We don't require all introns to match, just that any overlapping introns are compatible
    for t_intron in transcript_introns:
        for m_intron in model_introns:
            # Check if introns overlap
            if t_intron[0] <= m_intron[1] and t_intron[1] >= m_intron[0]:
                # If they overlap, they should be identical or very close
                # Allow small differences (e.g., 5 bp) to account for alignment uncertainty
                tolerance = 5
                if (
                    abs(t_intron[0] - m_intron[0]) > tolerance
                    or abs(t_intron[1] - m_intron[1]) > tolerance
                ):
                    return False

    return True


def select_best_utrs(utr_exons_list, strand, min_length=10, max_length=2000):
    """
    Select the best UTR exons from multiple transcript evidence.

    This function implements several strategies for selecting the most representative
    UTR structure from multiple transcript alignments.

    Parameters:
    -----------
    utr_exons_list : list
        List of lists, where each inner list contains UTR exon tuples (start, end)
    strand : str
        Strand of the gene model ('+' or '-')
    min_length : int, optional
        Minimum total length for a UTR to be considered (default: 10)
    max_length : int, optional
        Maximum total length for a UTR to be considered (default: 2000)

    Returns:
    --------
    tuple
        (best_utrs, method_used)
        - best_utrs: List of (start, end) tuples representing the best UTR exons
        - method_used: String describing the method used to select the UTRs
    """
    if not utr_exons_list:
        return None, "none"

    # Filter UTRs by length constraints
    filtered_utrs = []
    for utrs in utr_exons_list:
        total_length = sum(exon[1] - exon[0] + 1 for exon in utrs)
        if min_length <= total_length <= max_length:
            filtered_utrs.append(utrs)

    if not filtered_utrs:
        return None, "none"

    # If only one transcript has UTRs, use it
    if len(filtered_utrs) == 1:
        return sorted(filtered_utrs[0]), "single_transcript"

    # Calculate various metrics for each UTR set
    utr_metrics = []
    for i, utrs in enumerate(filtered_utrs):
        total_length = sum(exon[1] - exon[0] + 1 for exon in utrs)
        num_exons = len(utrs)
        utr_metrics.append(
            {
                "index": i,
                "utrs": utrs,
                "total_length": total_length,
                "num_exons": num_exons,
            }
        )

    # STRATEGY 1: Median length UTR
    # Sort by total length and take the median
    sorted_by_length = sorted(utr_metrics, key=lambda x: x["total_length"])
    median_idx = len(sorted_by_length) // 2
    median_utrs = sorted_by_length[median_idx]["utrs"]

    # STRATEGY 2: Most common exon count
    # Find the most common number of exons
    exon_counts = {}
    for metric in utr_metrics:
        count = metric["num_exons"]
        if count not in exon_counts:
            exon_counts[count] = 1
        else:
            exon_counts[count] += 1

    most_common_count = max(exon_counts.items(), key=lambda x: x[1])[0]

    # Filter UTRs with the most common exon count
    common_exon_utrs = [m["utrs"] for m in utr_metrics if m["num_exons"] == most_common_count]

    # If multiple UTRs have the most common exon count, take the median length one
    if len(common_exon_utrs) > 1:
        common_lengths = [sum(exon[1] - exon[0] + 1 for exon in utrs) for utrs in common_exon_utrs]
        median_length_idx = len(common_lengths) // 2
        common_exon_median_utrs = common_exon_utrs[median_length_idx]
    else:
        common_exon_median_utrs = common_exon_utrs[0]

    # STRATEGY 3: Consensus UTR structure
    # This is more complex - we try to find consensus exon boundaries
    # First, identify all unique exon boundaries
    all_starts = []
    all_ends = []
    for utrs in filtered_utrs:
        for exon in utrs:
            all_starts.append(exon[0])
            all_ends.append(exon[1])

    # Find consensus boundaries (those that appear in multiple transcripts)
    start_counts = {}
    end_counts = {}
    for start in all_starts:
        if start not in start_counts:
            start_counts[start] = 1
        else:
            start_counts[start] += 1

    for end in all_ends:
        if end not in end_counts:
            end_counts[end] = 1
        else:
            end_counts[end] += 1

    # Consider a boundary consensus if it appears in at least 1/3 of transcripts
    min_count = max(1, len(filtered_utrs) // 3)
    consensus_starts = [start for start, count in start_counts.items() if count >= min_count]
    consensus_ends = [end for end, count in end_counts.items() if count >= min_count]

    # Create consensus exons from consensus boundaries
    consensus_exons = []
    if strand == "+":
        # Sort boundaries
        consensus_starts.sort()
        consensus_ends.sort()

        # Match starts with ends to form exons
        for i, start in enumerate(consensus_starts):
            # Find the closest end that is after this start
            valid_ends = [end for end in consensus_ends if end > start]
            if valid_ends:
                closest_end = min(valid_ends)
                consensus_exons.append((start, closest_end))
    else:  # strand == '-'
        # For minus strand, we want to match from the end backwards
        consensus_starts.sort(reverse=True)
        consensus_ends.sort(reverse=True)

        # Match ends with starts to form exons
        for i, end in enumerate(consensus_ends):
            # Find the closest start that is before this end
            valid_starts = [start for start in consensus_starts if start < end]
            if valid_starts:
                closest_start = max(valid_starts)
                consensus_exons.append((closest_start, end))

    # Check if consensus approach yielded valid exons
    if consensus_exons:
        # Sort and remove overlaps
        consensus_exons.sort(key=lambda x: x[0])
        non_overlapping = [consensus_exons[0]]
        for exon in consensus_exons[1:]:
            if exon[0] > non_overlapping[-1][1]:
                non_overlapping.append(exon)

        # Check length constraints
        total_consensus_length = sum(exon[1] - exon[0] + 1 for exon in non_overlapping)
        if min_length <= total_consensus_length <= max_length:
            return non_overlapping, "consensus"

    # Choose the best strategy based on available data
    # If we have a good consensus, use that
    if consensus_exons and len(non_overlapping) > 0:
        return sorted(non_overlapping), "consensus"

    # If most transcripts agree on exon count, use that with median length
    if max(exon_counts.values()) > len(filtered_utrs) // 2:
        return sorted(common_exon_median_utrs), "common_exon_count"

    # Otherwise, use the median length UTR
    return sorted(median_utrs), "median_length"


def extend_utrs(
    consensus_models,
    transcripts,
    genome_fasta,
    min_utr_length=10,
    max_utr_length=2000,
    log=sys.stderr.write,
):
    """
    Extend consensus gene models with UTRs based on transcript evidence.

    This function examines transcript alignments that match consensus gene models
    and extends the models with 5' and 3' UTRs if supported by the evidence.
    Only transcripts with compatible intron/exon boundaries are considered.
    Properly handles spliced UTRs (UTRs containing introns).

    UTR extension is limited to avoid overlapping with neighboring genes on the same strand.

    Parameters:
    -----------
    consensus_models : dict
        Dictionary of consensus gene models
    transcripts : dict
        Dictionary of transcript alignments
    genome_fasta : str
        Path to the genome FASTA file
    min_utr_length : int, optional
        Minimum length for a UTR to be added (default: 10)
    max_utr_length : int, optional
        Maximum length for a UTR extension (default: 2000)
    log : callable, optional
        Function to use for logging (default: sys.stderr.write)

    Returns:
    --------
    dict
        Dictionary of consensus gene models with UTRs added where supported
    """

    # Track statistics
    models_with_5utr = 0
    models_with_3utr = 0
    total_transcripts = 0
    compatible_transcripts = 0
    utr_methods_used = defaultdict(int)

    # Build gene boundary index to prevent overlapping UTR extensions
    gene_boundaries = defaultdict(
        lambda: defaultdict(list)
    )  # {contig: {strand: [(start, end, gene_id)]}}
    for gene_id, model in consensus_models.items():
        if "location" in model and "strand" in model and "contig" in model:
            contig = model["contig"]
            strand = model["strand"]
            start, end = model["location"]
            gene_boundaries[contig][strand].append((start, end, gene_id))

    # Sort gene boundaries by position for efficient neighbor lookup
    for contig in gene_boundaries:
        for strand in gene_boundaries[contig]:
            gene_boundaries[contig][strand].sort(key=lambda x: x[0])

    def get_neighbor_boundaries(contig, strand, gene_id, current_start, current_end):
        """Get the boundaries of neighboring genes to limit UTR extension."""
        if contig not in gene_boundaries or strand not in gene_boundaries[contig]:
            return None, None

        genes_on_strand = gene_boundaries[contig][strand]
        upstream_boundary = None
        downstream_boundary = None

        for start, end, gid in genes_on_strand:
            if gid == gene_id:
                continue  # Skip the current gene

            # For plus strand: upstream is before gene start, downstream is after gene end
            # For minus strand: upstream is after gene end, downstream is before gene start
            if strand == "+":
                if end < current_start:  # Gene is upstream
                    upstream_boundary = max(upstream_boundary or 0, end)
                elif start > current_end:  # Gene is downstream
                    downstream_boundary = min(downstream_boundary or float("inf"), start)
            else:  # minus strand
                if start > current_end:  # Gene is upstream (in transcriptional direction)
                    upstream_boundary = min(upstream_boundary or float("inf"), start)
                elif end < current_start:  # Gene is downstream (in transcriptional direction)
                    downstream_boundary = max(downstream_boundary or 0, end)

        return upstream_boundary, downstream_boundary

    # Build separate InterLap objects for plus and minus strand transcripts by contig
    plus_transcript_interlap = defaultdict(interlap.InterLap)
    minus_transcript_interlap = defaultdict(interlap.InterLap)

    # Organize transcripts into InterLap objects by contig and strand
    for t_id, t in transcripts.items():
        contig = t["contig"]
        strand = t["strand"]
        t_min = min([x[0] for x in t["coords"]])
        t_max = max([x[1] for x in t["coords"]])

        # Store all needed information directly in the InterLap object
        if strand == "+":
            plus_transcript_interlap[contig].add((t_min, t_max, t_id, t["coords"]))
        else:
            minus_transcript_interlap[contig].add((t_min, t_max, t_id, t["coords"]))

    # Process each consensus model
    for model_id, model in list(consensus_models.items()):
        if "coords" not in model:
            continue

        # Get model coordinates and strand
        coords = model["coords"]
        strand = model["strand"]
        contig = model["contig"]

        # First, copy the original CDS coordinates
        model["cds"] = coords.copy()

        # Check if the model already has UTRs - if so, skip UTR extension
        if "five_prime_utr" in model or "three_prime_utr" in model:
            continue

        # Skip if no transcripts for this contig
        if contig not in plus_transcript_interlap and contig not in minus_transcript_interlap:
            continue

        # Get gene model boundaries
        g_min = min([x[0] for x in coords])
        g_max = max([x[1] for x in coords])

        # Use the appropriate InterLap object based on strand
        transcript_interlap = (
            plus_transcript_interlap[contig] if strand == "+" else minus_transcript_interlap[contig]
        )

        # Use InterLap to find overlapping transcripts on the same strand
        overlapping_transcripts = list(transcript_interlap.find((g_min, g_max)))
        total_transcripts += len(overlapping_transcripts)

        # Find transcript alignments that are compatible with this model
        matching_transcripts = []

        for _, _, t_id, t_coords in overlapping_transcripts:
            # Check if transcript has compatible intron/exon boundaries
            if check_intron_compatibility(coords, t_coords, strand):
                compatible_transcripts += 1

                # Calculate overlap score to ensure good match
                overlap_score = score_evidence(coords, t_coords)
                if overlap_score > 5:  # Only use transcripts with good overlap
                    matching_transcripts.append((t_id, t_coords))

        if not matching_transcripts:
            continue

        # Sort coordinates for consistent processing
        sorted_coords = sorted(coords, key=lambda x: x[0])

        # Get CDS boundaries (first and last exon)
        model_start = sorted_coords[0][0]
        model_end = sorted_coords[-1][1]

        # Initialize UTR lists
        five_prime_utrs = []
        three_prime_utrs = []

        # Process each matching transcript to find UTR exons
        for _, t_coords in matching_transcripts:
            t_coords_sorted = sorted(t_coords, key=lambda x: x[0])

            if strand == "+":
                # Find 5' UTR exons (before first CDS exon)
                five_prime_exons = []
                for exon in t_coords_sorted:
                    # Exon is entirely before the first CDS exon
                    if exon[1] < model_start:
                        # Check distance limit - don't extend 5' UTR too far upstream
                        # Also check for large gaps that suggest this is a different gene
                        gap_to_cds = model_start - exon[1]
                        if gap_to_cds <= max_utr_length and gap_to_cds <= 500:  # Max 500bp gap
                            five_prime_exons.append(exon)
                    # Exon overlaps with first CDS exon - partial UTR
                    elif exon[0] < model_start and exon[1] >= model_start:
                        # Check distance limit for the UTR portion
                        if model_start - exon[0] <= max_utr_length:
                            five_prime_exons.append((exon[0], model_start - 1))

                # Find 3' UTR exons (after last CDS exon)
                three_prime_exons = []
                for exon in t_coords_sorted:
                    # Exon is entirely after the last CDS exon
                    if exon[0] > model_end:
                        # Check distance limit - don't extend 3' UTR too far downstream
                        # Also check for large gaps that suggest this is a different gene
                        gap_to_cds = exon[0] - model_end
                        if gap_to_cds <= max_utr_length and gap_to_cds <= 500:  # Max 500bp gap
                            three_prime_exons.append(exon)
                    # Exon overlaps with last CDS exon - partial UTR
                    elif exon[0] <= model_end and exon[1] > model_end:
                        # Check distance limit for the UTR portion
                        if exon[1] - model_end <= max_utr_length:
                            three_prime_exons.append((model_end + 1, exon[1]))

                # Add to UTR collections if found
                if five_prime_exons:
                    five_prime_utrs.append(five_prime_exons)
                if three_prime_exons:
                    three_prime_utrs.append(three_prime_exons)
            else:
                # For minus strand, the logic is reversed
                # Find 3' UTR exons (before first CDS exon in genomic coordinates)
                three_prime_exons = []
                for exon in t_coords_sorted:
                    # Exon is entirely before the first CDS exon
                    if exon[1] < model_start:
                        # Check distance limit - don't extend 3' UTR too far upstream
                        # Also check for large gaps that suggest this is a different gene
                        gap_to_cds = model_start - exon[1]
                        if gap_to_cds <= max_utr_length and gap_to_cds <= 500:  # Max 500bp gap
                            three_prime_exons.append(exon)
                    # Exon overlaps with first CDS exon - partial UTR
                    elif exon[0] < model_start and exon[1] >= model_start:
                        # Check distance limit for the UTR portion
                        if model_start - exon[0] <= max_utr_length:
                            three_prime_exons.append((exon[0], model_start - 1))

                # Find 5' UTR exons (after last CDS exon in genomic coordinates)
                five_prime_exons = []
                for exon in t_coords_sorted:
                    # Exon is entirely after the last CDS exon
                    if exon[0] > model_end:
                        # Check distance limit - don't extend 5' UTR too far downstream
                        # Also check for large gaps that suggest this is a different gene
                        gap_to_cds = exon[0] - model_end
                        if gap_to_cds <= max_utr_length and gap_to_cds <= 500:  # Max 500bp gap
                            five_prime_exons.append(exon)
                    # Exon overlaps with last CDS exon - partial UTR
                    elif exon[0] <= model_end and exon[1] > model_end:
                        # Check distance limit for the UTR portion
                        if exon[1] - model_end <= max_utr_length:
                            five_prime_exons.append((model_end + 1, exon[1]))

                # Add to UTR collections if found
                if five_prime_exons:
                    five_prime_utrs.append(five_prime_exons)
                if three_prime_exons:
                    three_prime_utrs.append(three_prime_exons)

        # Flag to track if UTRs were added
        has_utrs = False

        # Get current gene model boundaries
        gene_start, gene_end = model["location"]

        # Create a copy of the original coords to modify
        merged_coords = coords.copy()

        # Process 5' UTRs if any exist
        if five_prime_utrs:
            # Select the best UTR using our new function
            best_five_prime_utrs, method_used = select_best_utrs(
                five_prime_utrs,
                strand,
                min_length=min_utr_length,
                max_length=max_utr_length,
            )

            utr_methods_used[f"5prime_{method_used}"] += 1

            if best_five_prime_utrs:
                has_utrs = True

                # Store the 5' UTR exons as a list of tuples - KEEP ORIGINAL COORDINATES
                model["five_prime_utr"] = best_five_prime_utrs

                # Get neighbor boundaries to limit UTR extension
                # Use original gene coordinates for neighbor detection
                orig_start, orig_end = model["location"]
                upstream_boundary, downstream_boundary = get_neighbor_boundaries(
                    contig, strand, model_id, orig_start, orig_end
                )

                # Update gene boundaries with neighbor checking
                if strand == "+":
                    proposed_start = min(exon[0] for exon in best_five_prime_utrs)
                    # For plus strand, 5' UTR extends upstream (lower coordinates)
                    if upstream_boundary is not None:
                        # Don't extend past the upstream neighbor (leave 2 bp gap minimum)
                        proposed_start = max(proposed_start, upstream_boundary + 2)
                    gene_start = min(gene_start, proposed_start)
                else:
                    proposed_end = max(exon[1] for exon in best_five_prime_utrs)
                    # For minus strand, 5' UTR extends upstream (higher coordinates)
                    if upstream_boundary is not None:
                        # Don't extend past the upstream neighbor (leave 2 bp gap minimum)
                        proposed_end = min(proposed_end, upstream_boundary - 2)
                    gene_end = max(gene_end, proposed_end)

                # Merge UTRs with coords for mRNA structure
                if strand == "+":
                    # For plus strand, check if first UTR exon connects to first CDS exon
                    sorted_utrs = sorted(best_five_prime_utrs, key=lambda x: x[0])
                    sorted_cds = sorted(coords, key=lambda x: x[0])

                    # Check if the last UTR exon connects to the first CDS exon
                    if sorted_utrs[-1][1] + 1 == sorted_cds[0][0]:
                        # Find the index of the first CDS exon in merged_coords
                        # Use a safer approach to find the index
                        first_cds_exon = sorted_cds[0]
                        found = False
                        for i, exon in enumerate(merged_coords):
                            if exon[0] == first_cds_exon[0] and exon[1] == first_cds_exon[1]:
                                # For merged coords, extend the first CDS exon to include the UTR
                                merged_coords[i] = (
                                    sorted_utrs[-1][0],
                                    first_cds_exon[1],
                                )
                                found = True
                                break

                        # Add all other UTR exons except the last one
                        if found:
                            for utr_exon in sorted_utrs[:-1]:
                                merged_coords.append(utr_exon)
                        else:
                            # If we couldn't find the exact exon, add all UTR exons
                            for utr_exon in sorted_utrs:
                                merged_coords.append(utr_exon)
                    else:
                        # Add all UTR exons as separate exons
                        for utr_exon in sorted_utrs:
                            merged_coords.append(utr_exon)
                else:
                    # For minus strand, check if first UTR exon connects to last CDS exon
                    sorted_utrs = sorted(best_five_prime_utrs, key=lambda x: x[0])
                    sorted_cds = sorted(coords, key=lambda x: x[0])

                    # Check if the first UTR exon connects to the last CDS exon
                    if sorted_cds[-1][1] + 1 == sorted_utrs[0][0]:
                        # Find the index of the last CDS exon in merged_coords
                        last_cds_exon = sorted_cds[-1]
                        found = False
                        for i, exon in enumerate(merged_coords):
                            if exon[0] == last_cds_exon[0] and exon[1] == last_cds_exon[1]:
                                # For merged coords, extend the last CDS exon to include the UTR
                                merged_coords[i] = (last_cds_exon[0], sorted_utrs[0][1])
                                found = True
                                break

                        # Add all other UTR exons except the first one
                        if found:
                            for utr_exon in sorted_utrs[1:]:
                                merged_coords.append(utr_exon)
                        else:
                            # If we couldn't find the exact exon, add all UTR exons
                            for utr_exon in sorted_utrs:
                                merged_coords.append(utr_exon)
                    else:
                        # Add all UTR exons as separate exons
                        for utr_exon in sorted_utrs:
                            merged_coords.append(utr_exon)

                models_with_5utr += 1

        # Process 3' UTRs if any exist
        if three_prime_utrs:
            # Select the best UTR using our new function
            best_three_prime_utrs, method_used = select_best_utrs(
                three_prime_utrs,
                strand,
                min_length=min_utr_length,
                max_length=max_utr_length,
            )

            utr_methods_used[f"3prime_{method_used}"] += 1

            if best_three_prime_utrs:
                has_utrs = True

                # Store the 3' UTR exons as a list of tuples - KEEP ORIGINAL COORDINATES
                model["three_prime_utr"] = best_three_prime_utrs

                # Get neighbor boundaries to limit UTR extension
                # Use original gene coordinates for neighbor detection
                orig_start, orig_end = model["location"]
                upstream_boundary, downstream_boundary = get_neighbor_boundaries(
                    contig, strand, model_id, orig_start, orig_end
                )

                # Update gene boundaries with neighbor checking
                if strand == "+":
                    proposed_end = max(exon[1] for exon in best_three_prime_utrs)
                    # For plus strand, 3' UTR extends downstream (higher coordinates)
                    if downstream_boundary is not None:
                        # Don't extend past the downstream neighbor (leave 2 bp gap minimum)
                        proposed_end = min(proposed_end, downstream_boundary - 2)
                    gene_end = max(gene_end, proposed_end)
                else:
                    proposed_start = min(exon[0] for exon in best_three_prime_utrs)
                    # For minus strand, 3' UTR extends downstream (lower coordinates)
                    if downstream_boundary is not None:
                        # Don't extend past the downstream neighbor (leave 2 bp gap minimum)
                        proposed_start = max(proposed_start, downstream_boundary + 2)
                    gene_start = min(gene_start, proposed_start)

                # Merge UTRs with coords for mRNA structure
                if strand == "+":
                    # For plus strand, check if first UTR exon connects to last CDS exon
                    sorted_utrs = sorted(best_three_prime_utrs, key=lambda x: x[0])
                    sorted_cds = sorted(coords, key=lambda x: x[0])

                    # Check if the first UTR exon connects to the last CDS exon
                    if sorted_cds[-1][1] + 1 == sorted_utrs[0][0]:
                        # Find the index of the last CDS exon in merged_coords
                        last_cds_exon = sorted_cds[-1]
                        found = False
                        for i, exon in enumerate(merged_coords):
                            if exon[0] == last_cds_exon[0] and exon[1] == last_cds_exon[1]:
                                # For merged coords, extend the last CDS exon to include the UTR
                                merged_coords[i] = (last_cds_exon[0], sorted_utrs[0][1])
                                found = True
                                break

                        # Add all other UTR exons except the first one
                        if found:
                            for utr_exon in sorted_utrs[1:]:
                                merged_coords.append(utr_exon)
                        else:
                            # If we couldn't find the exact exon, add all UTR exons
                            for utr_exon in sorted_utrs:
                                merged_coords.append(utr_exon)
                    else:
                        # Add all UTR exons as separate exons
                        for utr_exon in sorted_utrs:
                            merged_coords.append(utr_exon)
                else:
                    # For minus strand, check if last UTR exon connects to first CDS exon
                    sorted_utrs = sorted(best_three_prime_utrs, key=lambda x: x[0])
                    sorted_cds = sorted(coords, key=lambda x: x[0])

                    # Check if the last UTR exon connects to the first CDS exon
                    if sorted_utrs[-1][1] + 1 == sorted_cds[0][0]:
                        # Find the index of the first CDS exon in merged_coords
                        first_cds_exon = sorted_cds[0]
                        found = False
                        for i, exon in enumerate(merged_coords):
                            if exon[0] == first_cds_exon[0] and exon[1] == first_cds_exon[1]:
                                # For merged coords, extend the first CDS exon to include the UTR
                                merged_coords[i] = (
                                    sorted_utrs[-1][0],
                                    first_cds_exon[1],
                                )
                                found = True
                                break

                        # Add all other UTR exons except the last one
                        if found:
                            for utr_exon in sorted_utrs[:-1]:
                                merged_coords.append(utr_exon)
                        else:
                            # If we couldn't find the exact exon, add all UTR exons
                            for utr_exon in sorted_utrs:
                                merged_coords.append(utr_exon)
                    else:
                        # Add all UTR exons as separate exons
                        for utr_exon in sorted_utrs:
                            merged_coords.append(utr_exon)

                models_with_3utr += 1

        # Update gene location to include all UTRs
        if has_utrs:
            # Calculate final gene coordinates based on all actual features (UTRs + CDS)
            all_coords = []

            # Add CDS coordinates
            all_coords.extend(coords)

            # Add UTR coordinates if they exist
            if "five_prime_utr" in model:
                all_coords.extend(model["five_prime_utr"])
            if "three_prime_utr" in model:
                all_coords.extend(model["three_prime_utr"])

            # Calculate final gene boundaries from all actual coordinates
            final_start = min(coord[0] for coord in all_coords)
            final_end = max(coord[1] for coord in all_coords)

            model["location"] = (final_start, final_end)
            model["has_utrs"] = True
            # Update the coords with the merged version
            model["coords"] = sorted(merged_coords, key=lambda x: x[0])

        # Update model in the dictionary
        consensus_models[model_id] = model

        # Update the gene boundary index with the new coordinates after UTR extension
        if has_utrs and "location" in model:
            contig = model["contig"]
            strand = model["strand"]
            new_start, new_end = model["location"]

            # Find and update this gene's entry in the boundary index
            if contig in gene_boundaries and strand in gene_boundaries[contig]:
                for i, (start, end, gid) in enumerate(gene_boundaries[contig][strand]):
                    if gid == model_id:
                        gene_boundaries[contig][strand][i] = (new_start, new_end, gid)
                        break
                # Re-sort after updating
                gene_boundaries[contig][strand].sort(key=lambda x: x[0])

    log(
        f"Found {total_transcripts} overlapping transcripts, {compatible_transcripts} with compatible intron/exon boundaries"
    )
    log(
        f"Added 5' UTRs to {models_with_5utr} gene models and 3' UTRs to {models_with_3utr} gene models"
    )

    # Log the methods used for UTR selection
    if utr_methods_used:
        utr_methods = {}
        for method, count in sorted(utr_methods_used.items(), key=lambda x: x[1], reverse=True):
            if "none" in method:
                continue
            utr_methods[method] = count
        log(f"UTR selection methods used:\n{json.dumps(utr_methods, indent=2)}")

    return consensus_models


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


def cluster_interlap_original(obj):
    """
    Original cluster_interlap function (before fix) for debugging.
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


def cluster_interlap(obj):
    """
    Cluster genomic features using the interlap.reduce function.

    This function takes an interlap object containing genomic features and clusters them
    based on their coordinates. Features that overlap or are adjacent to each other are
    grouped into the same cluster. The function then assigns the original feature data
    to each cluster, ensuring each gene model is only assigned to one locus.

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

    # Fix for interlap.reduce() bug: merge bins that share genes
    # This handles cases where genes span bin boundaries
    merged_bins = []
    assigned_genes = set()

    for b in bins:
        hits = list(obj.find(b))
        current_genes = set(hit[2] for hit in hits)  # Gene IDs in this bin

        # Check if any genes in this bin are already assigned to a previous bin
        overlapping_genes = current_genes & assigned_genes

        if overlapping_genes:
            # Merge with the last bin that contains overlapping genes
            if merged_bins:
                last_bin = merged_bins[-1]
                # Expand the locus boundaries to include this bin
                new_start = min(last_bin["locus"][0], b[0])
                new_end = max(last_bin["locus"][1], b[1])
                last_bin["locus"] = (new_start, new_end)

                # Add new genes to the existing bin
                for hit in hits:
                    gene_id = hit[2]
                    if gene_id not in assigned_genes:
                        last_bin["genes"].append(hit[2:])
                        assigned_genes.add(gene_id)
            else:
                # This shouldn't happen, but handle it gracefully
                unique_hits = []
                for hit in hits:
                    gene_id = hit[2]
                    if gene_id not in assigned_genes:
                        unique_hits.append(hit)
                        assigned_genes.add(gene_id)

                if unique_hits:
                    merged_bins.append(
                        {
                            "locus": b,
                            "genes": [x[2:] for x in unique_hits],
                            "proteins": [],
                            "transcripts": [],
                            "repeats": [],
                        }
                    )
        else:
            # No overlapping genes, create a new bin
            unique_hits = []
            for hit in hits:
                gene_id = hit[2]
                if gene_id not in assigned_genes:
                    unique_hits.append(hit)
                    assigned_genes.add(gene_id)

            if unique_hits:
                merged_bins.append(
                    {
                        "locus": b,
                        "genes": [x[2:] for x in unique_hits],
                        "proteins": [],
                        "transcripts": [],
                        "repeats": [],
                    }
                )

    return merged_bins


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
                    info["location"],  # Add full gene location (includes UTRs)
                    info.get("5UTR", []),  # Add 5' UTR information
                    info.get("3UTR", []),  # Add 3' UTR information
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
                    info["location"],  # Add full gene location (includes UTRs)
                    info.get("5UTR", []),  # Add 5' UTR information
                    info.get("3UTR", []),  # Add 3' UTR information
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
        # Handle different gene tuple formats
        if len(gene) == 4:
            name, source, coords, cstart = gene
        elif len(gene) == 5:
            name, source, coords, cstart, full_location = gene
        elif len(gene) == 7:
            name, source, coords, cstart, full_location, utr_5, utr_3 = gene
        elif len(gene) == 9:
            start, end, name, source, coords, cstart, full_location, utr_5, utr_3 = gene
        else:
            raise ValueError(f"Unexpected gene tuple length: {len(gene)}")
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
            # Handle different gene tuple formats
            if len(gene) == 4:
                name, source, coords, cstart = gene
            elif len(gene) == 5:
                name, source, coords, cstart, full_location = gene
            elif len(gene) == 7:
                # Old format with 7 elements
                name, source, coords, cstart, full_location, utr_5, utr_3 = gene
            elif len(gene) == 9:
                # New format with 9 elements (includes UTR information)
                (
                    start,
                    end,
                    name,
                    source,
                    coords,
                    cstart,
                    full_location,
                    utr_5,
                    utr_3,
                ) = gene
            else:
                raise ValueError(f"Unexpected gene tuple length: {len(gene)}")
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
        # Handle different gene tuple formats
        first_gene = locus["genes"][0]
        if len(first_gene) == 4:
            name, source, coords, cstart = first_gene
            results = {name: {"source": source, "coords": coords, "score": 0}}
        elif len(first_gene) == 5:
            name, source, coords, cstart, full_location = first_gene
            results = {name: {"source": source, "coords": coords, "score": 0}}
        elif len(first_gene) == 7:
            name, source, coords, cstart, full_location, utr_5, utr_3 = first_gene
            results = {name: {"source": source, "coords": coords, "score": 0}}
        elif len(first_gene) == 9:
            start, end, name, source, coords, cstart, full_location, utr_5, utr_3 = first_gene
            results = {name: {"source": source, "coords": coords, "score": 0}}
        else:
            raise ValueError(f"Unexpected gene tuple length: {len(first_gene)}")
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
    # Handle both old format (4 elements) and new format (5 elements)
    sources = []
    for x in locus["genes"]:
        if len(x) >= 4 and x[1] not in derived:
            sources.append(x[1])  # source is at index 1
    if len(sources) > 0:
        src, n = Counter(sources).most_common(1)[0]
        if n > 1:
            # [['contig_76-snap.4', 'snap', [[(11710, 11713), (12068, 12155), (12543, 12881), (15021, 15918)]], 1]]
            seeds = []
            queries = []
            locs = []
            for x in locus["genes"]:
                # Handle both old format (4 elements) and new format (5 elements)
                if len(x) >= 4 and x[1] == src:
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
                    # Handle both old format (4 elements) and new format (5 elements)
                    if len(y) >= 4:
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
            # Handle different gene tuple formats
            if len(gene) == 4:
                name, source, coords, cstart = gene
                full_location = None
                utr_5 = []
                utr_3 = []
            elif len(gene) == 5:
                name, source, coords, cstart, full_location = gene
                utr_5 = []
                utr_3 = []
            elif len(gene) == 7:
                # Old format with 7 elements
                name, source, coords, cstart, full_location, utr_5, utr_3 = gene
            elif len(gene) == 9:
                # New format with 9 elements (includes UTR information)
                (
                    start,
                    end,
                    name,
                    source,
                    coords,
                    cstart,
                    full_location,
                    utr_5,
                    utr_3,
                ) = gene
            else:
                raise ValueError(f"Unexpected gene tuple length: {len(gene)}")
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
            "5UTR": utr_5,  # Add 5' UTR information
            "3UTR": utr_3,  # Add 3' UTR information
        }
        # Store full gene location (including UTRs) if available
        if full_location is not None:
            d["full_location"] = full_location
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
                        # Handle different gene tuple formats
                        if len(gene) == 4:
                            source = gene[1]  # Extract source
                        elif len(gene) == 5:
                            source = gene[1]  # Extract source
                        elif len(gene) == 7:
                            source = gene[1]  # Extract source
                        elif len(gene) == 9:
                            source = gene[
                                3
                            ]  # Extract source (different position in 9-element tuple)
                        else:
                            raise ValueError(f"Unexpected gene tuple length: {len(gene)}")
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
        # Handle different gene tuple formats to get the name
        if len(y) == 4:
            name = y[0]
        elif len(y) == 5:
            name = y[0]
        elif len(y) == 7:
            name = y[0]
        elif len(y) == 9:
            name = y[2]  # Different position in 9-element tuple
        else:
            raise ValueError(f"Unexpected gene tuple length: {len(y)}")
        dist[name] = {}
    for x in itertools.combinations(locus["genes"], 2):
        a, b = x
        # Handle different gene tuple formats
        if len(a) == 4 and len(b) == 4:
            a_name, a_coord = a[0], a[2]
            b_name, b_coord = b[0], b[2]
        elif len(a) == 5 and len(b) == 5:
            a_name, a_coord = a[0], a[2]
            b_name, b_coord = b[0], b[2]
        elif len(a) == 7 and len(b) == 7:
            a_name, a_coord = a[0], a[2]
            b_name, b_coord = b[0], b[2]
        elif len(a) == 9 and len(b) == 9:
            a_name, a_coord = a[2], a[4]  # Different positions in 9-element tuple
            b_name, b_coord = b[2], b[4]
        else:
            raise ValueError(f"Unexpected gene tuple lengths: {len(a)}, {len(b)}")
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
                sortedCDS = sorted(d["cds"], key=lambda tup: tup[0])
            else:
                sortedCoords = sorted(d["coords"], key=lambda tup: tup[0], reverse=True)
                sortedCDS = sorted(d["cds"], key=lambda tup: tup[0], reverse=True)
            # Check if location is valid
            try:
                g_start, g_end = d["location"]
            except (ValueError, TypeError):
                # Skip this gene if location is invalid
                continue

            # Preserve the original source name
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

            # if five_prime_utr then write here
            if "has_utrs" in d and d["has_utrs"] and "five_prime_utr" in d:
                for x in range(0, len(d["five_prime_utr"])):
                    futr_num = x + 1
                    exon = d["five_prime_utr"][x]
                    outfile.write(
                        "{}\tGFFtk\tfive_prime_UTR\t{}\t{}\t.\t{}\t.\tID={}.utr5p.{};Parent={}.mrna;\n".format(
                            d["contig"],
                            exon[0],
                            exon[1],
                            d["strand"],
                            locusTag,
                            futr_num,
                            locusTag,
                        )
                    )
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
            # if three_prime_utr then write here
            if "has_utrs" in d and d["has_utrs"] and "three_prime_utr" in d:
                for x in range(0, len(d["three_prime_utr"])):
                    t_num = x + 1
                    exon = d["three_prime_utr"][x]
                    outfile.write(
                        "{}\tGFFtk\tthree_prime_UTR\t{}\t{}\t.\t{}\t.\tID={}.utr5p.{};Parent={}.mrna;\n".format(
                            d["contig"],
                            exon[0],
                            exon[1],
                            d["strand"],
                            locusTag,
                            t_num,
                            locusTag,
                        )
                    )
            # now write the cds features
            if "cds" in d:
                for x in range(0, len(sortedCDS)):
                    outfile.write(
                        "{}\tGFFtk\tCDS\t{}\t{}\t.\t{}\t{}\tID={}.cds;Parent={}.mrna;\n".format(
                            d["contig"],
                            sortedCDS[x][0],
                            sortedCDS[x][1],
                            d["strand"],
                            current_phase,
                            locusTag,
                            locusTag,
                        )
                    )
                    current_phase = (
                        current_phase - (int(sortedCDS[x][1]) - int(sortedCDS[x][0]) + 1)
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

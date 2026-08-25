import concurrent.futures
import gzip
import io
import random
import re
import sys
import types
import uuid
from collections import OrderedDict
from urllib.parse import quote, unquote

from natsort import natsorted

from .fasta import RevComp, codon_table, fasta2dict, fasta2headers, getSeqRegions, translate
from .utils import zopen


def is_combined_gff_fasta(filename):
    """Check if a file contains both GFF3 and FASTA data.

    Parameters
    ----------
    filename : str
        Path to the file to check

    Returns
    -------
    bool
        True if file contains ##FASTA directive, False otherwise
    """
    if isinstance(filename, io.BytesIO):
        filename.seek(0)
        infile = filename
    else:
        infile = zopen(filename)

    try:
        for line in infile:
            if line.startswith("##FASTA"):
                return True
            # Stop checking after we see the first FASTA header or too many lines
            if line.startswith(">") or line.count("\n") > 10000:
                break
    finally:
        if not isinstance(filename, (io.BytesIO, io.StringIO)):
            infile.close()

    return False


def split_combined_gff_fasta(filename):
    """Split a combined GFF3+FASTA file into separate GFF3 and FASTA components.

    Parameters
    ----------
    filename : str
        Path to the combined file

    Returns
    -------
    tuple
        (gff_content, fasta_content) as file-like objects
    """
    gff_lines = []
    fasta_lines = []
    in_fasta_section = False

    if isinstance(filename, io.BytesIO):
        filename.seek(0)
        infile = filename
    else:
        infile = zopen(filename)

    try:
        for line in infile:
            if line.startswith("##FASTA"):
                in_fasta_section = True
                continue

            if in_fasta_section:
                fasta_lines.append(line)
            else:
                gff_lines.append(line)
    finally:
        if not isinstance(filename, (io.BytesIO, io.StringIO)):
            infile.close()

    # Create file-like objects from the content
    gff_content = io.StringIO("".join(gff_lines))
    fasta_content = io.StringIO("".join(fasta_lines))

    return gff_content, fasta_content


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


def _process_introns_to_exons(Genes, introns_by_parent, idParent, mrna_coordinates=None):
    """Process intron features to generate exon coordinates.

    For each mRNA that has introns, calculate exon coordinates by:
    1. Finding the mRNA boundaries (start, end)
    2. Sorting introns by position
    3. Generating exon coordinates from non-intronic regions

    Parameters
    ----------
    Genes : dict
        The genes dictionary being built
    introns_by_parent : dict
        Dictionary mapping parent IDs to list of intron coordinates
    idParent : dict
        Dictionary mapping feature IDs to their parent gene IDs
    mrna_coordinates : dict, optional
        Dictionary mapping transcript IDs to their mRNA feature coordinates
    """
    if mrna_coordinates is None:
        mrna_coordinates = {}
    for parent_id, intron_coords in introns_by_parent.items():
        # Find the gene that contains this mRNA
        gene_id = idParent.get(parent_id)
        if not gene_id or gene_id not in Genes:
            continue

        # Find the transcript index for this parent_id
        try:
            transcript_idx = Genes[gene_id]["ids"].index(parent_id)
            # Safety check for reasonable transcript index bounds
            if transcript_idx < 0 or transcript_idx > 100:  # Reasonable upper bound
                continue
        except (ValueError, KeyError):
            continue

        # Get the mRNA boundaries - we need to find the mRNA feature coordinates
        # For now, we'll calculate from the intron coordinates
        if not intron_coords:
            continue

        # Sort introns by start position
        intron_coords.sort(key=lambda x: x[0])

        # We need to find the actual mRNA boundaries
        # Priority order: mRNA features > existing exon coordinates > gene boundaries
        transcript_start = None
        transcript_end = None

        # First check if we have mRNA coordinates from mRNA features
        if parent_id in mrna_coordinates:
            # Use mRNA feature boundaries (highest priority)
            transcript_start, transcript_end = mrna_coordinates[parent_id]
        elif (
            transcript_idx < len(Genes[gene_id]["mRNA"]) and Genes[gene_id]["mRNA"][transcript_idx]
        ):
            # Use existing mRNA exon boundaries
            mrna_coords = Genes[gene_id]["mRNA"][transcript_idx]
            transcript_start = min(coord[0] for coord in mrna_coords)
            transcript_end = max(coord[1] for coord in mrna_coords)
        else:
            # Fallback to gene location as mRNA boundaries
            transcript_start, transcript_end = Genes[gene_id]["location"]

        # Safety check - ensure we have valid boundaries
        if transcript_start is None or transcript_end is None:
            continue

        # Generate exon coordinates by excluding intron regions
        exon_coords = []
        current_pos = transcript_start

        for intron_start, intron_end in intron_coords:
            # Add exon before this intron (if there's space)
            if current_pos < intron_start:
                exon_coords.append((current_pos, intron_start - 1))
            # Move past this intron
            current_pos = intron_end + 1

        # Add final exon after last intron (if there's space)
        if current_pos <= transcript_end:
            exon_coords.append((current_pos, transcript_end))

        # Store the calculated exon coordinates in the mRNA field
        # Ensure mRNA list exists and is properly sized
        if "mRNA" not in Genes[gene_id]:
            Genes[gene_id]["mRNA"] = []

        # Safely extend the mRNA list to accommodate this transcript index
        # Use direct extension instead of while loop to avoid infinite loops
        if transcript_idx >= 0 and transcript_idx < 100:  # Reasonable bounds check
            # Extend list directly to required size
            required_size = transcript_idx + 1
            current_size = len(Genes[gene_id]["mRNA"])
            if current_size < required_size:
                # Extend with empty lists
                Genes[gene_id]["mRNA"].extend([[] for _ in range(required_size - current_size)])

            # Set the mRNA coordinates with calculated exons
            Genes[gene_id]["mRNA"][transcript_idx] = exon_coords


def _handle_mrna_without_exons(Genes, mrna_coordinates):
    """Handle SGD-style mRNA features without explicit exon features.

    For genes where mRNA coordinates are empty (no exon features were found),
    but we have mRNA transcript information, use the stored mRNA feature boundaries
    as single-exon coordinates.

    This addresses SGD format where mRNA features exist but no explicit exon features.

    Parameters
    ----------
    Genes : dict
        The genes dictionary to process
    mrna_coordinates : dict
        Dictionary mapping transcript IDs to their mRNA feature coordinates
    """
    for gene_id, gene_data in Genes.items():
        # Check each transcript in this gene
        for i, transcript_coords in enumerate(gene_data.get("mRNA", [])):
            # If mRNA coordinates are empty but we have transcript IDs
            if not transcript_coords and i < len(gene_data.get("ids", [])):
                transcript_id = gene_data["ids"][i]

                # Use stored mRNA feature coordinates if available
                if transcript_id in mrna_coordinates:
                    mrna_start, mrna_end = mrna_coordinates[transcript_id]
                    gene_data["mRNA"][i] = [(mrna_start, mrna_end)]
                else:
                    # Fallback to gene location
                    gene_start, gene_end = gene_data.get("location", (0, 0))
                    if gene_start and gene_end:
                        gene_data["mRNA"][i] = [(gene_start, gene_end)]


def _subtract_introns_from_region(region, introns):
    """Subtract intron coordinates from a genomic region to get exonic parts.

    Parameters
    ----------
    region : tuple
        (start, end) coordinates of the region
    introns : list
        List of (start, end) intron coordinates

    Returns
    -------
    list
        List of (start, end) coordinates for exonic parts
    """
    if not introns:
        return [region]

    region_start, region_end = region
    exonic_parts = []

    # Sort introns by start position
    sorted_introns = sorted(introns, key=lambda x: x[0])

    current_pos = region_start

    for intron_start, intron_end in sorted_introns:
        # Skip introns that don't overlap with our region
        if intron_end < region_start or intron_start > region_end:
            continue

        # Adjust intron boundaries to region boundaries
        intron_start = max(intron_start, region_start)
        intron_end = min(intron_end, region_end)

        # Add exonic part before this intron
        if current_pos < intron_start:
            exonic_parts.append((current_pos, intron_start - 1))

        # Move past this intron
        current_pos = intron_end + 1

    # Add final exonic part after last intron
    if current_pos <= region_end:
        exonic_parts.append((current_pos, region_end))

    return exonic_parts


def _calculate_utrs_for_transcript(gene_data, transcript_idx, utr_introns=None):
    """Calculate and add UTR coordinates for a transcript.

    For SGD genes, UTRs are not explicitly annotated but can be calculated
    from the difference between mRNA and CDS coordinates. Accounts for introns
    within UTR regions.

    Parameters
    ----------
    gene_data : dict
        The gene data structure
    transcript_idx : int
        Index of the transcript to process
    utr_introns : dict, optional
        Dictionary with '5utr' and '3utr' keys containing intron coordinates
    """
    if utr_introns is None:
        utr_introns = {"5utr": [], "3utr": []}

    try:
        # Safety checks
        if (
            not gene_data
            or transcript_idx < 0
            or transcript_idx >= len(gene_data.get("mRNA", []))
            or transcript_idx >= len(gene_data.get("CDS", []))
        ):
            return

        mrna_coords = gene_data["mRNA"][transcript_idx]
        cds_coords = gene_data["CDS"][transcript_idx]

        if not mrna_coords or not cds_coords:
            return

        # Handle both single-exon and multi-exon transcripts
        if len(mrna_coords) >= 1 and len(cds_coords) >= 1:
            # Calculate overall mRNA boundaries from all exons
            mrna_start = min(coord[0] for coord in mrna_coords)
            mrna_end = max(coord[1] for coord in mrna_coords)

            # Find the overall CDS boundaries
            cds_start = min(coord[0] for coord in cds_coords)
            cds_end = max(coord[1] for coord in cds_coords)

            # Calculate UTRs
            utr5_coords = []
            utr3_coords = []

            # 5' UTR: from mRNA start to CDS start
            if mrna_start < cds_start:
                utr5_region = (mrna_start, cds_start - 1)
                # Subtract any 5' UTR introns to get exonic UTR coordinates
                utr5_coords = _subtract_introns_from_region(
                    utr5_region, utr_introns.get("5utr", [])
                )

            # 3' UTR: from CDS end to mRNA end
            if cds_end < mrna_end:
                utr3_region = (cds_end + 1, mrna_end)
                # Subtract any 3' UTR introns to get exonic UTR coordinates
                utr3_coords = _subtract_introns_from_region(
                    utr3_region, utr_introns.get("3utr", [])
                )

            # Ensure UTR lists exist and are properly sized
            if "5UTR" not in gene_data:
                gene_data["5UTR"] = []
            if "3UTR" not in gene_data:
                gene_data["3UTR"] = []

            # Extend lists safely without while loops
            max_needed = transcript_idx + 1

            # Extend 5UTR list if needed
            current_5utr_size = len(gene_data["5UTR"])
            if current_5utr_size < max_needed:
                gene_data["5UTR"].extend([[] for _ in range(max_needed - current_5utr_size)])

            # Extend 3UTR list if needed
            current_3utr_size = len(gene_data["3UTR"])
            if current_3utr_size < max_needed:
                gene_data["3UTR"].extend([[] for _ in range(max_needed - current_3utr_size)])

            # Store UTR coordinates
            gene_data["5UTR"][transcript_idx] = utr5_coords
            gene_data["3UTR"][transcript_idx] = utr3_coords

    except Exception:
        # Silently fail to avoid breaking the parser
        pass


def _process_deferred_cds(Genes, deferred_cds, idParent):
    """Process CDS features that were deferred due to feature order issues.

    This handles SGD files where CDS features appear before their parent mRNA features.

    Parameters
    ----------
    Genes : dict
        The genes dictionary being built
    deferred_cds : list
        List of CDS feature data that was deferred
    idParent : dict
        Dictionary mapping feature IDs to their parent gene IDs
    """
    for cds_data in deferred_cds:
        Parent = cds_data["Parent"]
        start = cds_data["start"]
        end = cds_data["end"]
        phase = cds_data["phase"]

        # Process CDS with multiple parents
        if "," in Parent:
            parents = Parent.split(",")
        else:
            parents = [Parent]

        for p in parents:
            if p in idParent:
                GeneFeature = idParent.get(p)
                if GeneFeature and GeneFeature in Genes:
                    # Find the transcript index for this parent
                    try:
                        i = Genes[GeneFeature]["ids"].index(p)
                        Genes[GeneFeature]["CDS"][i].append((start, end))
                        # Add phase
                        try:
                            phase_val = int(phase)
                        except (ValueError, TypeError):
                            phase_val = 0
                        Genes[GeneFeature]["phase"][i].append(phase_val)

                        # UTR calculation moved to after mRNA coordinate handling

                    except (ValueError, IndexError):
                        # Transcript not found, skip this CDS
                        continue


def _transfer_sgd_functional_annotations(gene_data, gene_attrs):
    """Transfer functional annotations from SGD gene features to transcript level.

    SGD stores rich functional annotations (GO terms, notes, dbxrefs) at the gene level
    but these should be propagated to transcript-level features for proper annotation.

    Parameters
    ----------
    gene_data : dict
        The gene data structure to update
    gene_attrs : dict
        Dictionary of gene-level attributes from SGD
    """
    # Extract functional annotations from gene attributes
    ontology_terms = gene_attrs.get("Ontology_term", "")
    note = gene_attrs.get("Note", "")
    dbxref = gene_attrs.get("dbxref", "")

    # Parse ontology terms (GO terms, SO terms, etc.)
    go_terms = []
    if ontology_terms:
        terms = ontology_terms.split(",")
        for term in terms:
            term = term.strip()
            if term.startswith("GO:"):
                go_terms.append(term)

    # Parse database cross-references
    db_xrefs = []
    if dbxref:
        # Handle multiple dbxrefs separated by commas
        xrefs = dbxref.split(",")
        for xref in xrefs:
            xref = xref.strip()
            if xref:
                db_xrefs.append(xref)

    # Parse and clean note field (URL decode)
    notes = []
    if note:
        # URL decode the note field
        import urllib.parse

        decoded_note = urllib.parse.unquote(note)
        # Replace %20 with spaces and %3B with semicolons, etc.
        decoded_note = decoded_note.replace("%20", " ").replace("%3B", ";").replace("%2C", ",")
        notes.append(decoded_note)

    # Transfer annotations to all transcripts in this gene
    num_transcripts = len(gene_data.get("type", []))
    for i in range(num_transcripts):
        # Ensure lists exist and are properly sized
        for field in ["go_terms", "db_xref", "note"]:
            if field not in gene_data:
                gene_data[field] = []

            # Extend list safely without while loop
            required_size = i + 1
            current_size = len(gene_data[field])
            if current_size < required_size:
                gene_data[field].extend([[] for _ in range(required_size - current_size)])

        # Add GO terms
        if go_terms:
            gene_data["go_terms"][i].extend(go_terms)

        # Add database cross-references
        if db_xrefs:
            gene_data["db_xref"][i].extend(db_xrefs)

        # Add notes
        if notes:
            gene_data["note"][i].extend(notes)


def _validate_sgd_gene_models(Genes, gene_attributes):
    """Validate and correct SGD gene models using so_term_name attribute.

    Uses the so_term_name attribute from gene features to validate and correct
    transcript types and ensure proper CDS assignment for protein-coding genes.

    Parameters
    ----------
    Genes : dict
        The genes dictionary to process
    gene_attributes : dict
        Dictionary mapping gene IDs to their parsed attributes
    """
    for gene_id, gene_data in Genes.items():
        if gene_id not in gene_attributes:
            continue

        attrs = gene_attributes[gene_id]
        so_term = attrs.get("so_term_name", "")

        # Handle protein-coding genes
        if so_term == "protein_coding_gene":
            # Correct transcript types from ncRNA to mRNA
            for i, transcript_type in enumerate(gene_data.get("type", [])):
                if transcript_type in ["ncRNA", "transcript"]:
                    gene_data["type"][i] = "mRNA"

            # Update products for mRNA transcripts
            gene_symbol = attrs.get("gene", gene_data.get("name", gene_id))
            for i, product in enumerate(gene_data.get("product", [])):
                if len(gene_data.get("type", [])) > i and gene_data["type"][i] == "mRNA":
                    if gene_symbol and gene_symbol != gene_id:
                        gene_data["product"][i] = f"{gene_symbol} protein"
                    else:
                        gene_data["product"][i] = "hypothetical protein"

            # Transfer functional annotations from gene to transcript level
            _transfer_sgd_functional_annotations(gene_data, attrs)


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
    # Collect introns to generate exon coordinates later
    introns_by_parent = {}  # {parent_id: [(start, end), ...]}
    # Collect UTR introns for proper UTR calculation
    utr_introns_by_parent = (
        {}
    )  # {parent_id: {'5utr': [(start, end), ...], '3utr': [(start, end), ...]}}
    # Collect mRNA feature coordinates for SGD-style parsing
    mrna_coordinates = {}  # {transcript_id: (start, end)}
    # Collect gene attributes for SGD validation
    gene_attributes = {}  # {gene_id: attributes_dict}
    # Collect CDS features for deferred processing (SGD feature order issue)
    deferred_cds = []  # [(line_data, info), ...]
    SeqRecords = fasta2headers(fasta)
    if isinstance(gff, (io.BytesIO, io.StringIO)):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    for line in infile:
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
            "snoRNA",
            "snRNA",
            "rRNA",
            "tmRNA",
            "pseudogene",
            "five_prime_UTR",
            "five_prime_utr",
            "three_prime_UTR",
            "three_prime_utr",
            "intron",
            "noncoding_exon",
            "pseudogenic_exon",
        ]:
            continue

        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Parent = None
        Name = None
        Product = None
        GeneFeature = None
        gbkey = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.split("=", 1)
                info[k] = v.strip()
            except (IndexError, ValueError):
                pass
        # now can lookup in info dict for values
        ID = info.get("ID", None)
        Parent = info.get("Parent", None)
        Name = info.get("Name", None)
        if "DBxref" in info:
            DBxref = info.get("DBxref", None)
        elif "Dbxref" in info:
            DBxref = info.get("Dbxref", None)
        elif "dbxref" in info:
            DBxref = info.get("dbxref", None)
        else:
            DBxref = None
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
            if attr not in [
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
                if attr not in errors["unparsed_attributes"]:
                    errors["unparsed_attributes"].append(attr)
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in [
            "gene",
            "pseudogene",
            "ncRNA_gene",
            "rRNA_gene",
            "snoRNA_gene",
            "snRNA_gene",
            "telomerase_RNA_gene",
            "transposable_element_gene",
            "tRNA_gene",
        ]:
            if ID not in Genes:
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
            # Store gene attributes for SGD validation
            gene_attributes[ID] = info
        else:
            if not Parent:
                errors["no_parent"].append(line)
                continue
            if feature in [
                "mRNA",
                "transcript",
                "tRNA",
                "ncRNA",
                "rRNA",
                "snRNA",
                "snoRNA",
                "tmRNA",
                "telomerase_RNA",
            ]:
                # required that we have an ID here, else it not parsable
                if not ID:
                    if Name:
                        ID = Name
                    else:
                        errors["no_id"].append(line)
                        continue
                if gbkey and gbkey == "misc_RNA":
                    feature = "ncRNA"
                if not Product:
                    if feature in ["mRNA", "transcript"]:
                        Product = "hypothetical protein"
                if Parent not in Genes:
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
                if ID not in idParent:
                    idParent[ID] = Parent
                # Store mRNA feature coordinates for SGD-style parsing
                mrna_coordinates[ID] = (start, end)
            # treat exon features (including non-standard exon types)
            elif feature in ["exon", "noncoding_exon", "pseudogenic_exon"]:
                if "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                # Defer CDS processing to handle SGD feature order issues
                # where CDS features appear before their parent mRNA features
                deferred_cds.append(
                    {
                        "contig": contig,
                        "source": source,
                        "feature": feature,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase,
                        "Parent": Parent,
                        "Name": Name,
                        "info": info,
                    }
                )
                continue
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
                        if GeneFeature not in Genes:
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
                        if GeneFeature not in Genes:
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
            # handle intron features - collect them for later processing
            elif feature == "intron":
                if Parent:
                    # Handle multiple parents (comma-separated)
                    parents = Parent.split(",") if "," in Parent else [Parent]
                    for p in parents:
                        p = p.strip()
                        if p not in introns_by_parent:
                            introns_by_parent[p] = []
                        introns_by_parent[p].append((start, end))
            # handle UTR intron features - collect them for UTR calculation
            elif feature in ["five_prime_UTR_intron", "three_prime_UTR_intron"]:
                if Parent:
                    # Handle multiple parents (comma-separated)
                    parents = Parent.split(",") if "," in Parent else [Parent]
                    for p in parents:
                        p = p.strip()
                        if p not in utr_introns_by_parent:
                            utr_introns_by_parent[p] = {"5utr": [], "3utr": []}

                        if feature == "five_prime_UTR_intron":
                            utr_introns_by_parent[p]["5utr"].append((start, end))
                        elif feature == "three_prime_UTR_intron":
                            utr_introns_by_parent[p]["3utr"].append((start, end))

    # Process introns to generate exon coordinates
    _process_introns_to_exons(Genes, introns_by_parent, idParent, mrna_coordinates)

    # Process deferred CDS features (SGD feature order issue)
    _process_deferred_cds(Genes, deferred_cds, idParent)

    # Post-processing: Handle SGD-style mRNA features without explicit exons
    # If mRNA coordinates are empty but we have mRNA features, use mRNA boundaries as exon coordinates
    _handle_mrna_without_exons(Genes, mrna_coordinates)

    # Calculate UTRs for SGD genes (after mRNA coordinates are populated)
    for gene_id, gene_data in Genes.items():
        for i in range(len(gene_data.get("mRNA", []))):
            # Get transcript ID for this index
            transcript_id = (
                gene_data.get("ids", [None])[i] if i < len(gene_data.get("ids", [])) else None
            )
            utr_introns = utr_introns_by_parent.get(transcript_id, {})
            _calculate_utrs_for_transcript(gene_data, i, utr_introns)

    # SGD-specific validation and correction using so_term_name
    _validate_sgd_gene_models(Genes, gene_attributes)

    if not isinstance(gff, (io.BytesIO, io.StringIO)):
        infile.close()
    return Genes, errors


def _gff_miniprot_parser(gff, fasta, Genes):
    # this is a specific parser for lh3 miniprot gff3 format; contains only mRNA, CDS, and stop_codon
    # CDS features do not have ID=, just Parent; stop_codon just has Parent and Rank
    # idea is to go through line by line and parse the records and add to Genes dictionary
    errors = {
        "contig_name": [],
        "columns": [],
        "comments": [],
        "unparsed_attributes": [],
        "no_parent": [],
        "no_id": [],
    }
    SeqRecords = fasta2headers(fasta)
    if isinstance(gff, io.BytesIO):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    for line in infile:
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
            "mRNA",
            "CDS",
            "stop_codon",
        ]:
            continue
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Name = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.split("=", 1)
                info[k] = v.strip()
            except (IndexError, ValueError):
                pass
        # now can lookup in info dict for values
        ID = info.get("ID", None)
        Parent = info.get("Parent", None)
        Target = info.get("Target", None)
        Identity = info.get("Identity", None)
        Positive = info.get("Positive", None)
        Rank = info.get("Rank", None)
        if Target:
            Name = Target.split()[0]
        # for error reporting capture unparsed keys
        for attr, value in info.items():
            if attr not in ["ID", "Parent", "Identity", "Rank", "Positive", "Target"]:
                if attr not in errors["unparsed_attributes"]:
                    errors["unparsed_attributes"].append(attr)
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in ["mRNA"]:
            if ID not in Genes:
                Genes[ID] = {
                    "name": Name,
                    "type": ["mRNA"],
                    "transcript": [],
                    "cds_transcript": [],
                    "protein": [],
                    "5UTR": [[]],
                    "3UTR": [[]],
                    "gene_synonym": [],
                    "codon_start": [],
                    "ids": [ID],
                    "CDS": [[]],
                    "mRNA": [[]],
                    "strand": strand,
                    "EC_number": [[]],
                    "location": (start, end),
                    "contig": contig,
                    "product": ["miniprot alignment"],
                    "source": source,
                    "phase": [[]],
                    "db_xref": [[]],
                    "go_terms": [[]],
                    "note": [
                        [
                            f"IDENTITY:{Identity}",
                            f"RANK:{Rank}",
                            f"Positive:{Positive}",
                        ]
                    ],
                    "partialStart": [],
                    "partialStop": [],
                    "pseudo": False,
                    "score": [],
                    "target": [],
                }
            else:
                if start < Genes[ID]["location"][0]:
                    Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                if end > Genes[ID]["location"][1]:
                    Genes[ID]["location"] = (Genes[ID]["location"][0], end)
                Genes[ID]["Note"].append(
                    [
                        f"IDENTITY:{Identity}",
                        f"RANK:{Rank}",
                        f"Positive:{Positive}",
                    ]
                )
        elif feature in ["CDS"]:
            if Parent:
                # determine which transcript this is get index from id
                i = Genes[Parent]["ids"].index(Parent)
                Genes[Parent]["CDS"][i].append((start, end))
                Genes[Parent]["mRNA"][i].append((start, end))
                Genes[Parent]["score"].append(round(float(Identity) * 100, 2))
                Genes[Parent]["target"].append(Target)
                # add phase
                try:
                    Genes[Parent]["phase"][i].append(int(phase))
                except ValueError:
                    Genes[Parent]["phase"][i].append("?")
        elif feature in ["stop_codon"]:
            if Parent:
                i = Genes[Parent]["ids"].index(Parent)
                # here we need to extend the last CDS if + and first if - strand
                if strand == "+":
                    lt = Genes[Parent]["CDS"][i][-1]
                    nLT = (lt[0], end)
                    Genes[Parent]["CDS"][i][-1] = nLT
                    Genes[Parent]["mRNA"][i][-1] = nLT
                elif strand == "-":
                    ft = Genes[Parent]["CDS"][i][0]
                    nFT = (start, ft[1])
                    Genes[Parent]["CDS"][i][0] = nFT
                    Genes[Parent]["mRNA"][i][0] = nFT
    if not isinstance(gff, (io.BytesIO, io.StringIO)):
        infile.close()
    return Genes, errors


def _gff_alignment_parser(gff, fasta, Genes):
    # this is a specific parser for EVM-like alignment evidence
    # features are effectively mRNA coords linked by identical ID
    # idea is to go through line by line and parse the records and add to Genes dictionary
    errors = {
        "contig_name": [],
        "columns": [],
        "comments": [],
        "unparsed_attributes": [],
        "no_parent": [],
        "no_id": [],
    }
    SeqRecords = fasta2headers(fasta)
    if isinstance(gff, io.BytesIO):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    for line in infile:
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
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Name = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.split("=", 1)
                info[k] = v.strip()
            except (IndexError, ValueError):
                pass
        # now can lookup in info dict for values
        ID = info.get("ID", None)
        Target = info.get("Target", None)
        if Target:
            Name = Target.split()[0]
        # for error reporting capture unparsed keys
        for attr, value in info.items():
            if attr not in ["ID", "Target"]:
                if attr not in errors["unparsed_attributes"]:
                    errors["unparsed_attributes"].append(attr)
        if feature in ["cDNA_match", "EST_match", "nucleotide_to_protein_match"]:
            if feature == "nucleotide_to_protein_match":
                product_name = "protein alignment"
            else:
                product_name = "transcript alignment"
            if ID not in Genes:
                Genes[ID] = {
                    "name": Name,
                    "type": ["ncRNA"],
                    "transcript": [],
                    "cds_transcript": [],
                    "protein": [],
                    "5UTR": [[]],
                    "3UTR": [[]],
                    "gene_synonym": [],
                    "codon_start": [],
                    "ids": [ID],
                    "CDS": [[]],
                    "mRNA": [[(start, end)]],
                    "strand": strand,
                    "EC_number": [[]],
                    "location": (start, end),
                    "contig": contig,
                    "product": [product_name],
                    "source": source,
                    "phase": [[]],
                    "db_xref": [[]],
                    "go_terms": [[]],
                    "note": [[]],
                    "partialStart": [],
                    "partialStop": [],
                    "pseudo": False,
                    "score": [round(float(score), 2)],
                    "target": [Target],
                }
            else:
                if start < Genes[ID]["location"][0]:
                    Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                if end > Genes[ID]["location"][1]:
                    Genes[ID]["location"] = (Genes[ID]["location"][0], end)
                Genes[ID]["mRNA"][0].append((start, end))
                Genes[ID]["score"].append(round(float(score), 2))
                Genes[ID]["target"].append(Target)
    if not isinstance(gff, (io.BytesIO, io.StringIO)):
        infile.close()
    return Genes, errors


# SGD parser temporarily removed due to implementation issues
def _gff_ncbi_parser(gff, fasta, Genes):
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
    # Collect introns to generate exon coordinates later
    introns_by_parent = {}  # {parent_id: [(start, end), ...]}
    SeqRecords = fasta2headers(fasta)
    if isinstance(gff, io.BytesIO):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    for line in infile:
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
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Parent = None
        Name = None
        Product = None
        GeneFeature = None
        gbkey = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.split("=", 1)
                info[k] = v.strip()
            except (IndexError, ValueError):
                pass
        # now can lookup in info dict for values
        ID = info.get("ID", None)
        if ID and ID.startswith(("gene-", "exon-", "rna-", "cds-")):
            ID = ID.split("-", 1)[1]
        Parent = info.get("Parent", None)
        if Parent and Parent.startswith(("gene-", "rna-")):
            Parent = Parent.split("-", 1)[1]
        Name = info.get("gene", None)
        if "DBxref" in info:
            DBxref = info.get("DBxref", None)
        elif "Dbxref" in info:
            DBxref = info.get("Dbxref", None)
        elif "dbxref" in info:
            DBxref = info.get("dbxref", None)
        else:
            DBxref = None
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
            if attr not in [
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
                if attr not in errors["unparsed_attributes"]:
                    errors["unparsed_attributes"].append(attr)
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in [
            "gene",
            "pseudogene",
            "ncRNA_gene",
            "rRNA_gene",
            "snoRNA_gene",
            "snRNA_gene",
            "telomerase_RNA_gene",
            "transposable_element_gene",
            "tRNA_gene",
        ]:
            if ID not in Genes:
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
            if feature in [
                "mRNA",
                "transcript",
                "tRNA",
                "ncRNA",
                "rRNA",
                "snoRNA",
                "snRNA",
                "telomerase_RNA",
            ]:
                if gbkey and gbkey == "misc_RNA":
                    feature = "ncRNA"
                if not Product:
                    if feature in ["mRNA", "transcript"]:
                        Product = "hypothetical protein"
                if Parent not in Genes:
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
                if ID not in idParent:
                    idParent[ID] = Parent
            # treat exon features (including non-standard exon types)
            elif feature in ["exon", "noncoding_exon", "pseudogenic_exon"]:
                if "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                        if GeneFeature not in Genes:
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
                                "db_xref": [DBxref],
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
                            if DBxref:
                                for dbx in DBxref:
                                    if dbx not in Genes[GeneFeature]["db_xref"][i]:
                                        Genes[GeneFeature]["db_xref"][i].append(dbx)
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
                        if GeneFeature not in Genes:
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
                        if GeneFeature not in Genes:
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
            # handle intron features - collect them for later processing
            elif feature == "intron":
                if Parent:
                    if Parent not in introns_by_parent:
                        introns_by_parent[Parent] = []
                    introns_by_parent[Parent].append((start, end))

    # Process introns to generate exon coordinates
    _process_introns_to_exons(Genes, introns_by_parent, idParent, {})

    if not isinstance(gff, (io.BytesIO, io.StringIO)):
        infile.close()
    return Genes, errors


def validate_models(annotation, fadict, logger=sys.stderr.write, table=1, gap_filter=False):
    if isinstance(logger, types.BuiltinFunctionType):
        log = logger
    else:
        log = logger.info
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
                    fadict,
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
            log(
                "ERROR in parsing gene {}\n{}\n{}\n".format(
                    gene, assert_lengths_fail, annotation[gene]
                )
            )
            raise SystemExit(1)
    return annotation


def validate_and_translate_models(
    k, v, SeqRecords, gap_filter=False, table=1, logger=sys.stderr.write
):
    if isinstance(logger, types.BuiltinFunctionType):
        log = logger
    else:
        log = logger.info
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
        "phase": [],
    }
    assert len(v["ids"]) == len(v["type"])
    for i in range(0, len(v["ids"])):
        if v["type"][i] in [
            "mRNA",
            "tRNA",
            "ncRNA",
            "rRNA",
            "transcript",
            "snoRNA",
            "snRNA",
            "tmRNA",
            "telomerase_RNA",
        ]:
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
                results["phase"].append([])
            else:
                # Sort CDS and phase together to maintain correspondence
                results["type"].append("mRNA")
                cds_phase_pairs = list(zip(v["CDS"][i], v["phase"][i]))
                if v["strand"] == "+":
                    sorted_pairs = sorted(cds_phase_pairs, key=lambda pair: pair[0][0])
                else:
                    sorted_pairs = sorted(
                        cds_phase_pairs, key=lambda pair: pair[0][0], reverse=True
                    )

                # Unzip the sorted pairs
                sortedCDS = [pair[0] for pair in sorted_pairs]
                sortedPhase = [pair[1] for pair in sorted_pairs]

                # get the codon_start from the first CDS phase + 1
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
                        numStops = protSeq.count("*")
                        if protSeq and protSeq[-1] == "*":
                            numStops -= 1
                        translateResults.append((y, numStops, protSeq))
                    sortedResults = sorted(translateResults, key=lambda tup: tup[1])
                    codon_start = sortedResults[0][0]
                    protSeq = sortedResults[0][2]
                    v["phase"][i] = codon_start - 1
                else:
                    try:
                        codon_start = int(sortedPhase[0]) + 1
                    except IndexError:
                        # Default to 1 if index is out of range
                        codon_start = 1
                    # translate and get protein sequence
                    protSeq = translate(cdsSeq, v["strand"], codon_start - 1, table=table)
                results["codon_start"].append(codon_start)
                if codon_start > 1:
                    if v["strand"] == "+":
                        cdsSeq = cdsSeq[codon_start - 1 :]
                    elif v["strand"] == "-":
                        endTrunc = len(cdsSeq) - codon_start - 1
                        cdsSeq = cdsSeq[0:endTrunc]
                    else:
                        log(f"ERROR nonsensical strand ({v['strand']}) for gene {v['ids'][i]}")
                results["cds_transcript"].append(cdsSeq)
                results["CDS"].append(sortedCDS)
                results["phase"].append(sortedPhase)
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
                    if protSeq.rstrip("*").count("*") > 0:
                        results["pseudo"] = True
        else:
            results["CDS"].append([])
            results["phase"].append([])
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


def _detect_format(gff):
    # this is incomplete search, but sniff if this is an NCBI GFF3 record
    parser = _gff_default_parser
    _format = "default"
    if isinstance(gff, (io.BytesIO, io.StringIO)):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    for line in infile:
        if line.startswith("#"):
            if "!processor NCBI annotwriter" in line:
                parser = _gff_ncbi_parser
                _format = "ncbi"
        else:
            break
    if not isinstance(gff, (io.BytesIO, io.StringIO)):
        infile.close()
    return parser, _format


def _detect_gtf_format(gff):
    # this is incomplete search, but sniff if this is an NCBI GFF3 record
    parser = _gtf_default_parser
    _format = "default"
    if isinstance(gff, io.BytesIO):
        gff.seek(0)
        infile = gff
    else:
        infile = zopen(gff)
    features = set()
    for i, line in enumerate(infile):
        if line.startswith(("\n", "#")):
            continue
        cols = line.split("\t")
        features.add(cols[2])
        if i == 1000:  # look at first 1000 lines should be sufficient?
            break
    if "exon" not in features:
        parser = _gtf_genemark_parser
        _format = "genemark"
    elif "gene" not in features:
        parser = _gtf_jgi_parser
        _format = "jgi"
    if not isinstance(gff, io.BytesIO):
        infile.close()
    return parser, _format


def _clean_ncbi_names(annot):
    Clean = {}
    for k, v in annot.items():
        new_ids = []
        for i, x in enumerate(v["ids"]):
            new_ids.append("{}-T{}".format(k, i + 1))
        v["ids"] = new_ids
        Clean[k] = v
    return Clean


def _longest_orf(annot, fadict, minlen=50, table=1):
    # this is for alignment phasing, where we go through each gene
    # and see if can find full length coding region
    Clean = {}
    for k, v in annot.items():
        if v["strand"] == "+":
            sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0])
        else:
            sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0], reverse=True)
        mrnaSeq = getSeqRegions(fadict, v["contig"], sortedExons)
        if v["strand"] == "-":
            searchSeq = RevComp(mrnaSeq).upper()
        else:
            searchSeq = mrnaSeq.upper()
        # get all possible ORFs from the mRNA sequence
        valid_starts = "|".join(codon_table[table]["start"])
        allORFS = re.findall(rf"((?:{valid_starts})(?:\S{{3}})*?T(?:AG|AA|GA))", searchSeq)
        longestORF = None
        # if you found ORFs, then we'll get the longest one and see if meets criterea
        if len(allORFS) > 0:
            longestORF = max(allORFS, key=len)
            if len(longestORF) > minlen * 3:  # than at least min prot length
                # prot = translate(longestORF, "+", 0, table=table)
                # now find out where this ORF is in the coords, sort all back to left --> right
                if v["strand"] == "-":
                    longestORF = RevComp(longestORF)
                    sortedExons = sorted(v["mRNA"][0], key=lambda tup: tup[0])
                # find left most position
                left_pos = mrnaSeq.upper().index(longestORF)
                # map the coords
                CDS = []
                cov = 0
                lenOrf = len(longestORF)
                c = 0
                for i, (s, e) in enumerate(sortedExons):
                    gap_len = e - s + 1
                    last = c
                    c += gap_len
                    if c > left_pos:
                        start_pos = s + left_pos - last
                        break

                for i, (s, f) in enumerate(sortedExons):
                    if s <= start_pos <= f:  # means should start in this exon
                        if (start_pos + lenOrf - 1) <= f:  # then single coord
                            CDS.append((start_pos, (start_pos + lenOrf - 1)))
                            break
                        else:  # spans multiple coords
                            CDS.append((start_pos, f))
                            cov += f - start_pos + 1
                    elif len(CDS) > 0:
                        if cov < lenOrf:
                            remainder = lenOrf - cov
                            if (f - s) < remainder:
                                CDS.append((s, f))
                                cov += f - s + 1
                            else:
                                CDS.append((s, s + remainder - 1))
                                break
                # see if this is correct
                cdsSeq = getSeqRegions(fadict, v["contig"], CDS)
                try:
                    assert cdsSeq.upper() == longestORF.upper()
                    # okay, then we can add CDS and change type
                    v["type"] = ["mRNA"]
                    if v["strand"] == "+":
                        v["CDS"] = [sorted(CDS, key=lambda tup: tup[0])]
                    else:
                        v["CDS"] = [sorted(CDS, key=lambda tup: tup[0], reverse=True)]
                    v["phase"] = [["?"] * len(CDS)]
                    v["partialStart"] = [False]
                    v["partialStop"] = [False]
                    Clean[k] = v
                except AssertionError:
                    v["partialStart"] = [False]
                    v["partialStop"] = [False]
                    Clean[k] = v
            else:
                v["partialStart"] = [False]
                v["partialStop"] = [False]
                Clean[k] = v
        else:  # did not find orfs so add as is
            v["partialStart"] = [False]
            v["partialStop"] = [False]
            Clean[k] = v
    return Clean


def _gtf_default_parser(gtf, fasta, Genes, gtf_format="default"):
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
    if isinstance(gtf, io.BytesIO):
        gtf.seek(0)
        infile = gtf
    else:
        infile = zopen(gtf)
    for line in infile:
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
            "transcript",
            "mRNA",
            "exon",
            "CDS",
            "three_prime_utr",
            "five_prime_utr",
            "5UTR",
            "3UTR",
        ]:
            continue
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Parent = None
        Name = None
        Product = None
        GeneFeature = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.rsplit(" ", 1)
                info[k.strip()] = v.strip().replace('"', "")
            except (IndexError, ValueError):
                pass
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in ["gene"]:
            ID = info.get("gene_id", None)
            Name = info.get("gene_name", None)
            feature_type = info.get("gene_biotype", None)
            if ID not in Genes:
                if feature_type and feature_type == "pseudogene":
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
                    "gene_synonym": [],
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
            if feature in ["mRNA", "transcript"]:
                ID = info.get("transcript_id", None)
                Parent = info.get("gene_id", None)
                Name = info.get("transcript_name", None)
                feature_type = info.get("transcript_biotype", None)
                if feature_type and feature_type == "protein_coding":
                    Product = "hypothetical protein"
                    display_feature = "mRNA"
                elif feature_type and feature_type in ["ncRNA", "tRNA"]:
                    display_feature = feature_type
                    Product = "hypothetical {}".format(feature_type)
                else:
                    display_feature = feature
                if Parent not in Genes:
                    Genes[Parent] = {
                        "name": Name,
                        "type": [display_feature],
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
                        "gene_synonym": [],
                        "db_xref": [],
                        "go_terms": [],
                        "EC_number": [],
                        "note": [],
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
                    Genes[Parent]["db_xref"].append([])
                    Genes[Parent]["EC_number"].append([])
                    Genes[Parent]["gene_synonym"] += []
                    Genes[Parent]["go_terms"].append([])
                    Genes[Parent]["note"].append([])
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
                if ID not in idParent:
                    idParent[ID] = Parent
            # treat exon features (including non-standard exon types)
            elif feature in ["exon", "noncoding_exon", "pseudogenic_exon"]:
                Parent = info.get("transcript_id", None)
                if Parent and "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    else:
                        GeneFeature - info.get("name", None)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                                "gene_synonym": [],
                            }
                        else:
                            # determine which transcript this is get index from id
                            i = Genes[GeneFeature]["ids"].index(p)
                            Genes[GeneFeature]["mRNA"][i].append((start, end))
            # treat codings sequence features
            elif feature == "CDS":
                ID = info.get("protein_id", None)
                Parent = info.get("transcript_id", None)
                if Parent and "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                                "phase": [[int(phase)]],
                                "db_xref": [],
                                "go_terms": [],
                                "EC_number": [],
                                "note": [],
                                "partialStart": [False],
                                "partialStop": [False],
                                "pseudo": False,
                                "gene_synonym": [],
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
            elif feature == "five_prime_utr" or feature == "5UTR":
                Parent = info.get("transcript_id", None)
                if Parent and "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                                "gene_synonym": [],
                            }
                        else:
                            # determine which transcript this is get index from id
                            i = Genes[GeneFeature]["ids"].index(p)
                            Genes[GeneFeature]["5UTR"][i].append((start, end))
            # treat 3' UTR
            elif feature == "three_prime_utr" or feature == "3UTR":
                Parent = info.get("transcript_id", None)
                if Parent and "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                                "gene_synonym": [],
                            }
                        else:
                            # determine which transcript this is get index from id
                            i = Genes[GeneFeature]["ids"].index(p)
                            Genes[GeneFeature]["3UTR"][i].append((start, end))
                # note we are ignore start_codon and stop_codon as its redudant and not needed
    if not isinstance(gtf, io.BytesIO):
        infile.close()
    return Genes, errors


def _gtf_genemark_parser(gtf, fasta, Genes, gtf_format="genemark"):
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
    if isinstance(gtf, io.BytesIO):
        gtf.seek(0)
        infile = gtf
    else:
        infile = zopen(gtf)
    for line in infile:
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
            "transcript",
            "mRNA",
            "exon",
            "CDS",
            "start_codon",
            "stop_codon",
            "three_prime_utr",
            "five_prime_utr",
            "5UTR",
            "3UTR",
        ]:
            continue
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        if source == "GeneMark.hmm3":
            source = "genemark"
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Parent = None
        Name = None
        Product = None
        GeneFeature = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.rsplit(" ", 1)
                info[k.strip()] = v.strip().replace('"', "")
            except (IndexError, ValueError):
                pass
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in ["gene"]:
            ID = info.get("gene_id", None)
            Name = info.get("gene_name", None)
            feature_type = info.get("gene_biotype", None)
            if ID not in Genes:
                if feature_type and feature_type == "pseudogene":
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
                    "gene_synonym": [],
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
            if feature in ["mRNA", "transcript"]:
                ID = info.get("transcript_id", None)
                Parent = info.get("gene_id", None)
                Name = info.get("transcript_name", None)
                feature_type = info.get("transcript_biotype", None)
                Product = "hypothetical protein"
                display_feature = "mRNA"
                if Parent not in Genes:
                    Genes[Parent] = {
                        "name": Name,
                        "type": [display_feature],
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
                        "gene_synonym": [],
                        "db_xref": [],
                        "go_terms": [],
                        "EC_number": [],
                        "note": [],
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
                    Genes[Parent]["db_xref"].append([])
                    Genes[Parent]["EC_number"].append([])
                    Genes[Parent]["gene_synonym"] += []
                    Genes[Parent]["go_terms"].append([])
                    Genes[Parent]["note"].append([])
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
                if ID not in idParent:
                    idParent[ID] = Parent
            # treat codings sequence features
            elif feature == "CDS":
                ID = info.get("protein_id", None)
                Parent = info.get("transcript_id", None)
                if "," in Parent:
                    parents = Parent.split(",")
                else:
                    parents = [Parent]
                for p in parents:
                    if p in idParent:
                        GeneFeature = idParent.get(p)
                    if GeneFeature:
                        if GeneFeature not in Genes:
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
                                "phase": [[int(phase)]],
                                "db_xref": [],
                                "go_terms": [],
                                "EC_number": [],
                                "note": [],
                                "partialStart": [False],
                                "partialStop": [False],
                                "pseudo": False,
                                "gene_synonym": [],
                            }
                        else:
                            # determine which transcript this is get index from id
                            i = Genes[GeneFeature]["ids"].index(p)
                            Genes[GeneFeature]["CDS"][i].append((start, end))
                            Genes[GeneFeature]["mRNA"][i].append((start, end))
                            # add phase
                            try:
                                Genes[GeneFeature]["phase"][i].append(int(phase))
                            except ValueError:
                                Genes[GeneFeature]["phase"][i].append("?")
    if not isinstance(gtf, io.BytesIO):
        infile.close()
    return Genes, errors


def _gtf_jgi_parser(gtf, fasta, Genes, gtf_format="jgi"):
    # this is a custom parser for old gff/gtf JGI format (stupid).
    # idea is to go through line by line and parse the records and add to Genes dictionary
    errors = {
        "contig_name": [],
        "columns": [],
        "comments": [],
        "unparsed_attributes": [],
        "no_parent": [],
        "no_id": [],
    }
    SeqRecords = fasta2headers(fasta)
    if isinstance(gtf, io.BytesIO):
        gtf.seek(0)
        infile = gtf
    else:
        infile = zopen(gtf)
    for line in infile:
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
            "exon",
            "CDS",
        ]:
            continue
        if contig not in SeqRecords:
            errors["contig_name"].append(line)
            continue
        attributes = unquote(attributes)
        source = unquote(source)
        feature = unquote(feature)
        start = int(start)
        end = int(end)
        ID = None
        Parent = None
        info = {}
        for field in attributes.split(";"):
            try:
                k, v = field.rsplit(" ", 1)
                info[k.strip()] = v.strip().replace('"', "")
            except (IndexError, ValueError):
                pass
        # we can get the ID
        ID = info.get("name", None)
        if ID is None:
            ID = info.get("gene_name", None)
        if ID is None:
            continue
        # now we can do add to dictionary these parsed values
        # genbank gff files are incorrect for tRNA so check if gbkey exists and make up gene on the fly
        if feature in ["exon"]:
            Parent = info.get("transcriptId", None)
            if Parent is None:
                Parent = info.get("transcript_id", None)
            if Parent is None:
                continue
            if ID not in Genes:
                Genes[ID] = {
                    "name": None,
                    "type": [],
                    "transcript": [],
                    "cds_transcript": [],
                    "protein": [],
                    "5UTR": [[]],
                    "3UTR": [[]],
                    "gene_synonym": [],
                    "codon_start": [],
                    "ids": [f"{ID}-T1"],
                    "CDS": [[]],
                    "mRNA": [[(start, end)]],
                    "strand": strand,
                    "EC_number": [[]],
                    "location": (start, end),
                    "contig": contig,
                    "product": [],
                    "source": source,
                    "phase": [[]],
                    "db_xref": [[]],
                    "go_terms": [[]],
                    "note": [[]],
                    "partialStart": [],
                    "partialStop": [],
                    "pseudo": False,
                }
            else:
                if start < Genes[ID]["location"][0]:
                    Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                if end > Genes[ID]["location"][1]:
                    Genes[ID]["location"] = (Genes[ID]["location"][0], end)
                Genes[ID]["mRNA"][0].append((start, end))
        else:  # then its a CDS feature
            Parent = info.get("proteinId", None)
            if Parent is None:
                Parent = info.get("protein_id", None)
            if Parent is None:
                continue
            if ID not in Genes:
                Genes[ID] = {
                    "name": None,
                    "type": ["mRNA"],
                    "transcript": [],
                    "cds_transcript": [],
                    "protein": [],
                    "5UTR": [[]],
                    "3UTR": [[]],
                    "gene_synonym": [],
                    "codon_start": [],
                    "ids": [f"{ID}-T1"],
                    "CDS": [[(start, end)]],
                    "mRNA": [[(start, end)]],
                    "strand": strand,
                    "EC_number": [[]],
                    "location": (start, end),
                    "contig": contig,
                    "product": ["hypothetical protein"],
                    "source": source,
                    "phase": [[int(phase)]],
                    "db_xref": [[]],
                    "go_terms": [[]],
                    "note": [[]],
                    "partialStart": [],
                    "partialStop": [],
                    "pseudo": False,
                }
            else:
                Genes[ID]["CDS"][0].append((start, end))
                if len(Genes[ID]["type"]) == 0:
                    Genes[ID]["type"].append("mRNA")
                if len(Genes[ID]["product"]) == 0:
                    Genes[ID]["product"].append("hypothetical protein")
                # add phase
                try:
                    Genes[ID]["phase"][0].append(int(phase))
                except ValueError:
                    Genes[ID]["phase"][0].append("?")

    if not isinstance(gtf, io.BytesIO):
        infile.close()
    return Genes, errors


def gtf2dict(
    gtf,
    fasta,
    annotation=False,
    table=1,
    debug=False,
    gap_filter=False,
    gtf_format="auto",
    logger=sys.stderr.write,
):
    """Convert GTF and FASTA to standardized GFFtk dictionary format.

    Annotation file in GTF format and genome FASTA file are parsed. The result is a dictionary that is keyed
    by locus_tag (gene name) and the value is a nested dictionary containing feature information.

    Parameters
    ----------
    gtf : filename : str
        annotation text file in GTF format
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
    # parser format
    if gtf_format == "auto":
        gtf_parser, _format = _detect_gtf_format(gtf)
        _format = "default"
    elif gtf_format == "genemark":
        gtf_parser = _gtf_genemark_parser
        _format = "genemark"
    annotation, parse_errors = gtf_parser(gtf, fasta, annotation, gtf_format=_format)
    SeqRecords = fasta2dict(fasta)
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
    if debug:
        logger(
            "Parsed {} contigs containing {} gene models\n".format(len(SeqRecords), len(annotation))
        )
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    annotation = validate_models(
        annotation, SeqRecords, logger=logger, table=table, gap_filter=gap_filter
    )

    return annotation


def gff2dict(
    gff,
    fasta,
    annotation=False,
    table=1,
    debug=False,
    gap_filter=False,
    gff_format="auto",
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

    # Check if this is a combined GFF3+FASTA file
    # This happens when fasta is None/False and the gff file contains both GFF3 and FASTA data
    if (fasta is None or fasta is False or gff == fasta) and is_combined_gff_fasta(gff):
        if debug:
            logger("Detected combined GFF3+FASTA file, splitting content\n")
        gff_content, fasta_content = split_combined_gff_fasta(gff)
        gff = gff_content
        fasta = fasta_content

    # autodetect format
    if gff_format == "auto":
        gff_parser, _format = _detect_format(gff)
    elif gff_format == "ncbi-euk":
        gff_parser = _gff_ncbi_parser
        _format = "ncbi"
    elif gff_format == "miniprot":
        gff_parser = _gff_miniprot_parser
        _format = "miniprot"
    elif gff_format == "alignment":
        gff_parser = _gff_alignment_parser
        _format = "alignment"
    # elif gff_format == "sgd":
    #     gff_parser = _gff_sgd_parser
    #     _format = "sgd"
    else:
        gff_parser = _gff_default_parser
        _format = "default"
    annotation, parse_errors = gff_parser(gff, fasta, annotation)
    SeqRecords = fasta2dict(fasta)
    if _format == "ncbi":  # clean up identifer names
        annotation = _clean_ncbi_names(annotation)
    if _format == "alignment":
        annotation = _longest_orf(annotation, SeqRecords)
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
                    print_errors = random.sample(err_v, min(10, len(err_v)))
                    logger("{}\n".format("\n".join(print_errors)))
    if debug:
        logger(
            "Parsed {} contigs containing {} gene models\n".format(len(SeqRecords), len(annotation))
        )
    # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
    annotation = validate_models(
        annotation, SeqRecords, logger=logger, table=table, gap_filter=gap_filter
    )

    return annotation


def simplifyGO(inputList):
    simple = []
    for x in inputList:
        if x.startswith("GO:"):
            simple.append(x.strip())
        elif " " in x:
            simple.append(x.split(" ")[1])
    return simple


def dict2gff3(infile, output=False, debug=False, source=False, newline=False, url_encode=False):
    """Convert GFFtk standardized annotation dictionary to GFF3 file.

    Annotation dictionary generated by gff2dict or tbl2dict passed as input. This function then write to GFF3 format

    Parameters
    ----------
    infile : dict of dict
        standardized annotation dictionary keyed by locus_tag
    output : str, default=sys.stdout
        annotation file in GFF3 format
    debug : bool, default=False
        print debug information to stderr
    source : str, default=False
        override source field in GFF3 output
    newline : bool, default=False
        add newline after each gene
    url_encode : bool, default=False
        URL encode attribute values for downstream tool compatibility

    """

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(infile.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    if output:
        if hasattr(output, "write"):  # It's a file-like object (StringIO, etc.)
            gffout = output
        elif output.endswith(".gz"):
            copen = gzip.open
            mopen = "wt"
            gffout = copen(output, mopen)
        else:
            copen = open
            mopen = "w"
            gffout = copen(output, mopen)
    else:
        gffout = sys.stdout
    gffout.write("##gff-version 3\n")
    for k, v in list(sortedGenes.items()):
        genefeature = "gene"
        if "pseudo" in v:
            if v["pseudo"]:
                genefeature = "pseudogene"
        # option to update/modify/change the source
        if source:
            new_source = source
        else:
            new_source = v["source"]
        if v["name"]:
            if "gene_synonym" in v and len(v["gene_synonym"]) > 0:
                gffout.write(
                    "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};Alias={:};\n".format(
                        v["contig"],
                        new_source,
                        genefeature,
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
                    "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Name={:};\n".format(
                        v["contig"],
                        new_source,
                        genefeature,
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
                    "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Alias={:};\n".format(
                        v["contig"],
                        new_source,
                        genefeature,
                        v["location"][0],
                        v["location"][1],
                        v["strand"],
                        k,
                        ",".join(v["gene_synonym"]),
                    )
                )
            else:
                gffout.write(
                    "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};\n".format(
                        v["contig"],
                        new_source,
                        genefeature,
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
                if v["type"][i] == "mRNA":
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0])
                if "5UTR" in v and v["5UTR"][i]:
                    sortedFive = sorted(v["5UTR"][i], key=lambda tup: tup[0])
                if "3UTR" in v and v["3UTR"][i]:
                    sortedThree = sorted(v["3UTR"][i], key=lambda tup: tup[0])
            else:
                sortedExons = sorted(v["mRNA"][i], key=lambda tup: tup[0], reverse=True)
                if v["type"][i] == "mRNA":
                    sortedCDS = sorted(v["CDS"][i], key=lambda tup: tup[0], reverse=True)
                if "5UTR" in v and v["5UTR"][i]:
                    sortedFive = sorted(v["5UTR"][i], key=lambda tup: tup[0], reverse=True)
                if "3UTR" in v and v["3UTR"][i]:
                    sortedThree = sorted(v["3UTR"][i], key=lambda tup: tup[0], reverse=True)
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
                        if "," not in values:
                            CleanedNote.append(base + ":" + values)
                        else:
                            for y in values.split(","):
                                CleanedNote.append(base + ":" + y)
                    else:
                        CleanedNote.append(x.replace(",", ""))
                extraAnnotations = extraAnnotations + "Note={:};".format(",".join(CleanedNote))
            # Calculate mRNA boundaries from exon coordinates for this transcript
            if v["mRNA"][i] and len(v["mRNA"][i]) > 0:
                # Use actual exon coordinates to determine mRNA boundaries
                try:
                    mrna_start = min(coord[0] for coord in v["mRNA"][i] if len(coord) >= 2)
                    mrna_end = max(coord[1] for coord in v["mRNA"][i] if len(coord) >= 2)
                except (ValueError, TypeError):
                    # Fallback if coordinate calculation fails
                    mrna_start = v["location"][0]
                    mrna_end = v["location"][1]
            else:
                # Fallback to gene location if no exons
                mrna_start = v["location"][0]
                mrna_end = v["location"][1]

            # Apply URL encoding if requested
            product_value = v["product"][i]
            extra_annotations = extraAnnotations
            if url_encode:
                product_value = quote(str(product_value), safe="")
                # URL encode values in extraAnnotations while preserving structure
                if extra_annotations:
                    # Split by semicolons, encode values after equals signs
                    encoded_parts = []
                    for part in extra_annotations.split(";"):
                        if "=" in part and part.strip():
                            key, value = part.split("=", 1)
                            encoded_value = quote(str(value), safe="")
                            encoded_parts.append(f"{key}={encoded_value}")
                        elif part.strip():
                            encoded_parts.append(part)
                    extra_annotations = ";".join(encoded_parts)

            # now write mRNA feature
            gffout.write(
                "{:}\t{:}\t{:}\t{:}\t{:}\t.\t{:}\t.\tID={:};Parent={:};product={:};{:}\n".format(
                    v["contig"],
                    new_source,
                    v["type"][i],
                    mrna_start,
                    mrna_end,
                    v["strand"],
                    v["ids"][i],
                    k,
                    product_value,
                    extra_annotations,
                )
            )
            if v["type"][i] in ["mRNA", "tRNA", "ncRNA", "tmRNA", "snRNA", "snoRNA"]:
                if "5UTR" in v and v["5UTR"][i]:
                    # if 5'UTR then write those first
                    num_5utrs = len(v["5UTR"][i])
                    if num_5utrs > 0:
                        for z in range(0, num_5utrs):
                            u_num = z + 1
                            gffout.write(
                                "{:}\t{:}\tfive_prime_UTR\t{:}\t{:}\t.\t{:}\t.\tID={:}.utr5p{:};Parent={:};\n".format(
                                    v["contig"],
                                    new_source,
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
                            new_source,
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
                                    new_source,
                                    sortedThree[z][0],
                                    sortedThree[z][1],
                                    v["strand"],
                                    v["ids"][i],
                                    u_num,
                                    v["ids"][i],
                                )
                            )
            if v["type"][i] == "mRNA":
                # if a pseudogene do not output CDS
                if "pseudo" in v:
                    if v["pseudo"]:
                        continue
                num_cds = len(v["CDS"][i])
                # GFF3 phase is 1 less than flat file
                current_phase = v["codon_start"][i] - 1
                for y in range(0, num_cds):
                    cds_num = y + 1
                    gffout.write(
                        "{:}\t{:}\tCDS\t{:}\t{:}\t.\t{:}\t{:}\tID={:}.cds{:};Parent={:};\n".format(
                            v["contig"],
                            new_source,
                            sortedCDS[y][0],
                            sortedCDS[y][1],
                            v["strand"],
                            current_phase,
                            v["ids"][i],
                            cds_num,
                            v["ids"][i],
                        )
                    )
                    current_phase = (
                        current_phase - (int(sortedCDS[y][1]) - int(sortedCDS[y][0]) + 1)
                    ) % 3
                    if current_phase == 3:
                        current_phase = 0
        if newline:
            gffout.write("\n")
    if output and not hasattr(
        output, "write"
    ):  # Only close if it's a filename, not a file-like object
        gffout.close()


def dict2gtf(infile, output=False, source=False):
    """Convert GFFtk standardized annotation dictionary to GTF file.

    Annotation dictionary generated by gff2dict or tbl2dict passed as input. This function
    then write to GTF format, notably this function only writes protein coding CDS features.

    Parameters
    ----------
    infile : dict of dict
        standardized annotation dictionary keyed by locus_tag
    output : str, default=sys.stdout
        annotation in GTF format

    """

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(infile.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    if output:
        if output.endswith(".gz"):
            copen = gzip.open
            mopen = "wt"
        else:
            copen = open
            mopen = "w"
        gtfout = copen(output, mopen)
    else:
        gtfout = sys.stdout
    for k, v in list(sortedGenes.items()):
        for i in range(0, len(v["ids"])):
            if len(v["CDS"][i]) < 1:
                continue
            # option to change source
            if source:
                new_source = source
            else:
                new_source = v["source"]
            # create attributes string
            attributes = 'gene_id "{:}"; transcript_id "{:}";'.format(k, v["ids"][i])
            if len(v["5UTR"][i]) > 0:
                for utr in v["5UTR"][i]:
                    gtfout.write(
                        "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                            v["contig"],
                            new_source,
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
                        new_source,
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
                            new_source,
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
                    current_phase = (current_phase - (int(cds[1]) - int(cds[0]) + 1)) % 3
                    if current_phase == 3:
                        current_phase = 0
                    if v["strand"] == "+":
                        if x == len(v["CDS"][i]) - 1:  # this is last one
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    new_source,
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
                                    new_source,
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
                                        new_source,
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
                                raise SystemExit(1)
                    else:
                        if x == len(v["CDS"][i]) - 1:  # this is last one
                            gtfout.write(
                                "{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n".format(
                                    v["contig"],
                                    new_source,
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
                                    new_source,
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
                                    new_source,
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
                            new_source,
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


def dict2combined_gff_fasta(annotation_dict, fasta_dict, output=False, debug=False, source=False):
    """Write GFFtk annotation dictionary and FASTA sequences to combined GFF3+FASTA format.

    Parameters
    ----------
    annotation_dict : dict
        GFFtk standardized annotation dictionary
    fasta_dict : dict
        Dictionary of sequences keyed by contig name
    output : str or file handle, default=False
        Output file path or handle. If False, writes to stdout
    debug : bool, default=False
        Print debug information
    source : str, default=False
        Override source field in GFF3 output

    Returns
    -------
    None
    """
    # First write the GFF3 part
    if output:
        if isinstance(output, str):
            outfile = open(output, "w")
        else:
            outfile = output
    else:
        outfile = sys.stdout

    # Write GFF3 header
    outfile.write("##gff-version 3\n")

    # Write sequence regions if we have fasta data
    if fasta_dict:
        for contig, seq in fasta_dict.items():
            outfile.write(f"##sequence-region {contig} 1 {len(seq)}\n")

    # Write the annotation data using existing dict2gff3 function
    # We'll capture the GFF3 output and write it
    import io

    gff_buffer = io.StringIO()
    dict2gff3(annotation_dict, output=gff_buffer, debug=debug, source=source)
    gff_content = gff_buffer.getvalue()
    gff_buffer.close()

    # Skip the ##gff-version 3 line since we already wrote it
    gff_lines = gff_content.split("\n")
    for line in gff_lines:
        if line.startswith("##gff-version"):
            continue
        if line.strip():  # Skip empty lines
            outfile.write(line + "\n")

    # Write the FASTA section
    if fasta_dict:
        outfile.write("##FASTA\n")
        for contig, seq in fasta_dict.items():
            outfile.write(f">{contig}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                outfile.write(seq[i : i + 80] + "\n")

    if output and isinstance(output, str):
        outfile.close()


def dict2gff3alignments(
    infile,
    output=False,
    debug=False,
    alignments="transcript",
    source=False,
    newline=False,
):
    """Convert GFFtk standardized annotation dictionary to GFF3 alignments file.

    Annotation dictionary generated by gff2dict or tbl2dict passed as input. Output format is GFF3-alignment, aka EVM evidence format

    Parameters
    ----------
    infile : dict of dict
        standardized annotation dictionary keyed by locus_tag
    output : str, default=sys.stdout
        annotation file in GFF3 format
    debug : bool, default=False
        print debug information to stderr

    """

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(infile.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    # then loop through and write GFF3 format
    if output:
        if output.endswith(".gz"):
            copen = gzip.open
            mopen = "wt"
        else:
            copen = open
            mopen = "w"
        gffout = copen(output, mopen)
    else:
        gffout = sys.stdout
    # set the feature type
    feature_type = "EST_match"
    if alignments == "protein":
        feature_type = "nucleotide_to_protein_match"
    gffout.write("##gff-version 3\n")
    # here we want to write each alignment as separate line
    for k, v in list(sortedGenes.items()):
        # print(k, v["mRNA"], v["score"], v["target"])
        if source:
            new_source = source
        else:
            new_source = v["source"]
        for i in range(0, len(v["mRNA"][0])):
            gffout.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\tID={};Target={};\n".format(
                    v["contig"],
                    new_source,
                    feature_type,
                    v["mRNA"][0][i][0],
                    v["mRNA"][0][i][1],
                    v["score"][i],
                    v["strand"],
                    k,
                    v["target"][i],
                )
            )
        if newline:
            gffout.write("\n")
    if output:
        gffout.close()

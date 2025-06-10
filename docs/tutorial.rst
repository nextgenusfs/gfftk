Tutorial
========

This tutorial will guide you through common tasks using gfftk.

Installation
-----------

You can install gfftk using pip:

.. code-block:: bash

    pip install gfftk

For more installation options, see the :doc:`installation guide <install>`.

Basic GFF3 Operations
--------------------

Parsing a GFF3 File
~~~~~~~~~~~~~~~~~~

Let's start by parsing a GFF3 file:

.. code-block:: python

    import gfftk

    # Parse a GFF3 file
    gff_dict = gfftk.gff.gff2dict("input.gff3", "genome.fasta")

    # Print the number of genes
    print(f"Number of genes: {len(gff_dict)}")

Modifying Gene Annotations
~~~~~~~~~~~~~~~~~~~~~~~~~

You can modify gene annotations in the parsed GFF3 data:

.. code-block:: python

    import gfftk

    # Parse a GFF3 file
    gff_dict = gfftk.gff.gff2dict("input.gff3", "genome.fasta")

    # Modify gene annotations
    for gene_id, gene in gff_dict.items():
        # Add a note to each gene
        if "note" not in gene:
            gene["note"] = []
        gene["note"].append("Modified by gfftk")

        # Update the source
        gene["source"] = "gfftk"

        # Update mRNA sources
        for mrna in gene.get("mRNA", []):
            mrna["source"] = "gfftk"

    # Write the modified data back to a GFF3 file
    gfftk.gff.dict2gff3(gff_dict, output="modified.gff3")

Filtering Genes
~~~~~~~~~~~~~

You can filter genes based on various criteria using both the Python API and command-line interface.

**Manual Filtering with Python API:**

.. code-block:: python

    import gfftk

    # Parse a GFF3 file
    gff_dict = gfftk.gff.gff2dict("input.gff3", "genome.fasta")

    # Filter genes by length
    filtered_genes = {}
    for gene_id, gene in gff_dict.items():
        gene_length = gene["location"][1] - gene["location"][0] + 1
        if gene_length >= 1000:  # Only keep genes >= 1000 bp
            filtered_genes[gene_id] = gene

    # Write the filtered data back to a GFF3 file
    gfftk.gff.dict2gff3(filtered_genes, output="filtered.gff3")

**Built-in Filtering with Convert Command:**

The convert command provides built-in filtering options using ``--grep`` and ``--grepv`` flags:

.. code-block:: bash

    # Keep only kinase genes
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.gff3 --grep product:kinase

    # Remove augustus predictions
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.gff3 --grepv source:augustus

    # Case-insensitive filtering
    gfftk convert -i input.gff3 -f genome.fasta -o results.gff3 --grep product:KINASE:i

    # Combined filtering: keep kinases but remove augustus predictions
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.gff3 \
        --grep product:kinase --grepv source:augustus

**Filter Pattern Syntax:**

- Basic pattern: ``key:pattern`` (e.g., ``product:kinase``)
- Case-insensitive: ``key:pattern:i`` (e.g., ``product:KINASE:i``)
- Regex patterns: ``key:regex_pattern`` (e.g., ``contig:^chr[0-9]+$``)
- Multiple patterns: Use multiple ``--grep`` or ``--grepv`` flags

**Common Filter Examples:**

.. code-block:: bash

    # Filter by gene product
    gfftk convert -i input.gff3 -f genome.fasta -o transporters.gff3 --grep product:transporter

    # Filter by annotation source
    gfftk convert -i input.gff3 -f genome.fasta -o genemark_only.gff3 --grep source:genemark

    # Filter by chromosome/contig
    gfftk convert -i input.gff3 -f genome.fasta -o chr1_genes.gff3 --grep contig:chr1

    # Filter by strand
    gfftk convert -i input.gff3 -f genome.fasta -o plus_strand.gff3 --grep strand:\\+

    # Remove hypothetical proteins
    gfftk convert -i input.gff3 -f genome.fasta -o known_proteins.gff3 \
        --grepv product:"hypothetical.*protein"

**Available Filter Keys:**

You can filter on any annotation attribute including:

- ``product`` - Gene product/function
- ``source`` - Annotation source (augustus, genemark, etc.)
- ``name`` - Gene name
- ``note`` - Gene notes/comments
- ``contig`` - Chromosome/contig name
- ``strand`` - DNA strand (+ or -)
- ``type`` - Feature type
- ``db_xref`` - Database cross-references
- ``go_terms`` - Gene Ontology terms

Format Conversion
---------------

Converting GFF3 to GTF
~~~~~~~~~~~~~~~~~~~~~

You can convert a GFF3 file to GTF format using the command line:

.. code-block:: bash

    gfftk convert -i input.gff3 -f genome.fasta -o output.gtf

Or using the Python API:

.. code-block:: python

    import gfftk

    # Convert GFF3 to GTF
    gfftk.convert.gff2gtf("input.gff3", "genome.fasta", "output.gtf")

**Converting with Filtering:**

You can combine format conversion with filtering:

.. code-block:: bash

    # Convert only kinase genes to GTF
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.gtf --grep product:kinase

    # Convert to GTF excluding augustus predictions
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.gtf --grepv source:augustus

Converting GFF3 to BED
~~~~~~~~~~~~~~~~~~~~~

You can convert a GFF3 file to BED format using the command line:

.. code-block:: bash

    gfftk convert -i input.gff3 -f bed -o output.bed

Or using the Python API:

.. code-block:: python

    import gfftk

    # Convert GFF3 to BED
    gfftk.convert.gff2bed("input.gff3", "output.bed")

Converting GFF3 to TBL
~~~~~~~~~~~~~~~~~~~~~

You can convert a GFF3 file to TBL format (for GenBank submission) using the command line:

.. code-block:: bash

    gfftk convert -i input.gff3 -f tbl -g genome.fasta -o output.tbl

Or using the Python API:

.. code-block:: python

    import gfftk

    # Convert GFF3 to TBL
    gfftk.convert.gff2tbl("input.gff3", "genome.fasta", "output.tbl")

Extracting Protein Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~

You can extract protein sequences from a GFF3 file using the command line:

.. code-block:: bash

    gfftk convert -i input.gff3 -f genome.fasta -o proteins.fasta --output-format proteins

Or using the Python API:

.. code-block:: python

    import gfftk

    # Extract protein sequences
    gfftk.convert.gff2proteins("input.gff3", "genome.fasta", "proteins.fasta")

**Extracting Filtered Protein Sequences:**

You can extract proteins for specific gene sets:

.. code-block:: bash

    # Extract only kinase proteins
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.faa \
        --output-format proteins --grep product:kinase

    # Extract proteins excluding hypothetical proteins
    gfftk convert -i input.gff3 -f genome.fasta -o known_proteins.faa \
        --output-format proteins --grepv product:"hypothetical.*protein"

Extracting Transcript Sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can extract transcript sequences from a GFF3 file using the command line:

.. code-block:: bash

    gfftk convert -i input.gff3 -f genome.fasta -o transcripts.fasta --output-format transcripts

Or using the Python API:

.. code-block:: python

    import gfftk

    # Extract transcript sequences
    gfftk.convert.gff2transcripts("input.gff3", "genome.fasta", "transcripts.fasta")

**Extracting Filtered Transcript Sequences:**

You can extract transcripts for specific gene sets:

.. code-block:: bash

    # Extract transcripts from specific chromosome
    gfftk convert -i input.gff3 -f genome.fasta -o chr1_transcripts.fasta \
        --output-format transcripts --grep contig:chr1

    # Extract transcripts from high-confidence predictions
    gfftk convert -i input.gff3 -f genome.fasta -o confident_transcripts.fasta \
        --output-format transcripts --grepv source:augustus

Consensus Gene Models
-------------------

Generating Consensus Gene Models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can generate consensus gene models from multiple sources using the command line:

.. code-block:: bash

    gfftk consensus -i input1.gff3 input2.gff3 -f genome.fasta -o consensus.gff3

Or using the Python API:

.. code-block:: python

    import gfftk

    # Generate consensus gene models
    consensus = gfftk.consensus.generate_consensus(
        ["input1.gff3", "input2.gff3"],
        "genome.fasta",
        weights={"input1": 1, "input2": 2},
        threshold=3,
    )

    # Write the consensus gene models to a GFF3 file
    gfftk.gff.dict2gff3(consensus, output="consensus.gff3")

Using Weights for Consensus Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can assign different weights to different input sources:

.. code-block:: bash

    gfftk consensus -i input1.gff3 input2.gff3 input3.gff3 -f genome.fasta -o consensus.gff3 -w weights.json

Where weights.json is a JSON file with the following structure:

.. code-block:: json

    {
        "input1": 1,
        "input2": 2,
        "input3": 3
    }

Advanced Topics
-------------

Working with GenBank Files
~~~~~~~~~~~~~~~~~~~~~~~~

You can convert between GFF3 and GenBank formats:

.. code-block:: python

    import gfftk

    # Convert GFF3 to TBL (for GenBank submission)
    gfftk.genbank.gff2tbl("input.gff3", "genome.fasta", "output.tbl")

    # Convert GFF3 to GenBank
    gfftk.genbank.gff2gbk("input.gff3", "genome.fasta", "output.gbk")

    # Convert GenBank to GFF3
    gfftk.genbank.gbk2gff("input.gbk", "output.gff3")

Comparing GFF3 Files
------------------

You can compare two GFF3 files to identify differences using the command line:

.. code-block:: bash

    gfftk compare -i input1.gff3 -c input2.gff3 -f genome.fasta -o comparison.txt

Or using the Python API:

.. code-block:: python

    import gfftk

    # Parse the GFF3 files
    gff_dict1 = gfftk.gff.gff2dict("input1.gff3", "genome.fasta")
    gff_dict2 = gfftk.gff.gff2dict("input2.gff3", "genome.fasta")

    # Compare the GFF3 files
    comparison = gfftk.compare.compareAnnotations(gff_dict1, gff_dict2, "genome.fasta")

    # Print the comparison results
    print(f"Shared genes: {len(comparison['shared'])}")
    print(f"Unique to input1: {len(comparison['unique1'])}")
    print(f"Unique to input2: {len(comparison['unique2'])}")

Working with FASTA Files
~~~~~~~~~~~~~~~~~~~~~

gfftk provides functions for working with FASTA files:

.. code-block:: python

    import gfftk

    # Parse a FASTA file
    fasta_dict = gfftk.fasta.fasta2dict("genome.fasta")

    # Get the length of each sequence
    for seq_id, seq in fasta_dict.items():
        print(f"{seq_id}: {len(seq)} bp")

    # Reverse complement a sequence
    rev_comp = gfftk.fasta.RevComp(fasta_dict["seq1"])

    # Translate a sequence
    protein = gfftk.fasta.translate(fasta_dict["seq1"], "+", 0)

    # Extract a region from a sequence
    region = gfftk.fasta.getSeqRegions(fasta_dict, [["seq1", 1, 100]])[0]

    # Write a FASTA file
    gfftk.fasta.dict2fasta(fasta_dict, "output.fasta")

GFF3 File Manipulation
~~~~~~~~~~~~~~~~~~~

gfftk provides several commands for manipulating GFF3 files:

1. **Sorting GFF3 Files**

.. code-block:: bash

    gfftk sort -i input.gff3 -o sorted.gff3

2. **Sanitizing GFF3 Files**

.. code-block:: bash

    gfftk sanitize -i input.gff3 -o sanitized.gff3

3. **Renaming Features in GFF3 Files**

.. code-block:: bash

    gfftk rename -i input.gff3 -o renamed.gff3 -p PREFIX

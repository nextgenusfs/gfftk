Usage
=====

Command-Line Interface
---------------------

gfftk provides a command-line interface with several subcommands for working with GFF3 files.

Basic Usage
~~~~~~~~~~

.. code-block:: bash

    gfftk [subcommand] [options]

Available Subcommands
~~~~~~~~~~~~~~~~~~~~

* ``consensus``: Generate consensus gene models from multiple sources
* ``convert``: Convert between file formats
* ``sort``: Sort GFF3 files
* ``sanitize``: Clean up GFF3 files
* ``rename``: Rename features in GFF3 files
* ``stats``: Generate statistics about GFF3 files
* ``compare``: Compare GFF3 files

Detailed Command Options
~~~~~~~~~~~~~~~~~~~~~

**consensus**

.. code-block:: bash

    gfftk consensus -h

    usage: gfftk consensus [-h] -i GFF3 [GFF3 ...] -f FASTA [-o OUTPUT] [-w WEIGHTS] [-t THRESHOLD] [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3 [GFF3 ...], --input GFF3 [GFF3 ...]
                            Input GFF3 file(s)
      -f FASTA, --fasta FASTA
                            Genome FASTA file
      -o OUTPUT, --output OUTPUT
                            Output GFF3 file
      -w WEIGHTS, --weights WEIGHTS
                            JSON file with weights for each input
      -t THRESHOLD, --threshold THRESHOLD
                            Score threshold for consensus
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**convert**

.. code-block:: bash

    gfftk convert -h

    usage: gfftk convert [-h] -i GFF3 -o OUTPUT -f {gtf,bed,tbl,proteins,transcripts} [-g GENOME] [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -o OUTPUT, --output OUTPUT
                            Output file
      -f {gtf,bed,tbl,proteins,transcripts}, --format {gtf,bed,tbl,proteins,transcripts}
                            Output format
      -g GENOME, --genome GENOME
                            Genome FASTA file (required for some formats)
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**sort**

.. code-block:: bash

    gfftk sort -h

    usage: gfftk sort [-h] -i GFF3 -o OUTPUT [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -o OUTPUT, --output OUTPUT
                            Output GFF3 file
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**sanitize**

.. code-block:: bash

    gfftk sanitize -h

    usage: gfftk sanitize [-h] -i GFF3 -o OUTPUT [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -o OUTPUT, --output OUTPUT
                            Output GFF3 file
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**rename**

.. code-block:: bash

    gfftk rename -h

    usage: gfftk rename [-h] -i GFF3 -o OUTPUT -p PREFIX [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -o OUTPUT, --output OUTPUT
                            Output GFF3 file
      -p PREFIX, --prefix PREFIX
                            Prefix for gene IDs
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**stats**

.. code-block:: bash

    gfftk stats -h

    usage: gfftk stats [-h] -i GFF3 [-o OUTPUT] [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -o OUTPUT, --output OUTPUT
                            Output file for statistics
      --no-progress        Disable progress bar
      --debug              Enable debug logging

**compare**

.. code-block:: bash

    gfftk compare -h

    usage: gfftk compare [-h] -i GFF3 -c GFF3 -f FASTA [-o OUTPUT] [--no-progress] [--debug]

    options:
      -h, --help            show this help message and exit
      -i GFF3, --input GFF3
                            Input GFF3 file
      -c GFF3, --compare GFF3
                            GFF3 file to compare against
      -f FASTA, --fasta FASTA
                            Genome FASTA file
      -o OUTPUT, --output OUTPUT
                            Output file for comparison results
      --no-progress        Disable progress bar
      --debug              Enable debug logging

Examples
~~~~~~~~

1. Generate consensus gene models:

.. code-block:: bash

    gfftk consensus -i input1.gff3 input2.gff3 -f genome.fasta -o consensus.gff3

2. Convert a GFF3 file to GTF format:

.. code-block:: bash

    gfftk convert -i input.gff3 -o output.gtf -f gtf

3. Sort a GFF3 file:

.. code-block:: bash

    gfftk sort -i input.gff3 -o sorted.gff3

4. Generate statistics about a GFF3 file:

.. code-block:: bash

    gfftk stats -i input.gff3

5. Compare two GFF3 files:

.. code-block:: bash

    gfftk compare -i input1.gff3 -c input2.gff3 -o comparison.txt

Python API
---------

gfftk can also be used as a Python library. The library provides a comprehensive set of functions for working with GFF3 files.

Core Modules
~~~~~~~~~~~

* ``gfftk.gff``: Functions for parsing, manipulating, and writing GFF3 files
* ``gfftk.consensus``: Functions for generating consensus gene models
* ``gfftk.convert``: Functions for converting between file formats
* ``gfftk.compare``: Functions for comparing GFF3 files
* ``gfftk.fasta``: Functions for working with FASTA files
* ``gfftk.genbank``: Functions for working with GenBank files

Basic Usage
~~~~~~~~~~

.. code-block:: python

    import gfftk

    # Parse a GFF3 file
    gff_dict = gfftk.gff.gff2dict("input.gff3", "genome.fasta")

    # Modify the GFF3 data
    # ...

    # Write the modified data back to a GFF3 file
    gfftk.gff.dict2gff3(gff_dict, output="output.gff3")

GFF3 Data Structure
~~~~~~~~~~~~~~~~~

The GFF3 data structure used by gfftk is a nested dictionary with the following structure:

.. code-block:: python

    {
        "gene_id": {
            "type": "gene",
            "contig": "contig_name",
            "location": [start, end],
            "strand": "+" or "-",
            "source": "source_name",
            "score": score_value,
            "phase": phase_value,
            "ID": "gene_id",
            "Name": "gene_name",
            # Other attributes...
            "mRNA": [
                {
                    "type": "mRNA",
                    "id": "mrna_id",
                    "parent": "gene_id",
                    "exon": [[start1, end1], [start2, end2], ...],
                    "CDS": [[start1, end1], [start2, end2], ...],
                    # Other features and attributes...
                },
                # More mRNAs...
            ],
        },
        # More genes...
    }

Examples
~~~~~~~~

1. Parse a GFF3 file and extract gene information:

.. code-block:: python

    import gfftk

    # Parse a GFF3 file
    gff_dict = gfftk.gff.gff2dict("input.gff3", "genome.fasta")

    # Print information about each gene
    for gene_id, gene in gff_dict.items():
        print(f"Gene ID: {gene_id}")
        print(f"Location: {gene['contig']}:{gene['location'][0]}-{gene['location'][1]}")
        print(f"Strand: {gene['strand']}")
        print(f"Number of mRNAs: {len(gene.get('mRNA', []))}")
        print()

2. Convert a GFF3 file to BED format:

.. code-block:: python

    import gfftk

    # Convert GFF3 to BED
    gfftk.convert.gff2bed("input.gff3", "output.bed")

3. Generate consensus gene models:

.. code-block:: python

    import gfftk

    # Generate consensus gene models
    consensus = gfftk.consensus.consensus_gene_models(
        ["input1.gff3", "input2.gff3"],
        "genome.fasta",
        weights={"input1": 1, "input2": 2},
        threshold=3,
    )

    # Write the consensus gene models to a GFF3 file
    gfftk.gff.dict2gff3(consensus, output="consensus.gff3")

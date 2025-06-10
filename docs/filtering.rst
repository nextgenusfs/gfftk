Filtering Annotations
====================

GFFtk provides powerful filtering capabilities through the ``--grep`` and ``--grepv`` options in the convert command. These options allow you to filter gene annotations based on any attribute using flexible pattern matching.

Basic Usage
-----------

Filter Syntax
~~~~~~~~~~~~~

The basic syntax for filtering is:

- ``--grep key:pattern`` - Keep annotations matching the pattern
- ``--grepv key:pattern`` - Remove annotations matching the pattern

Pattern Format
~~~~~~~~~~~~~~

Patterns can be specified in several formats:

- ``key:pattern`` - Basic string matching
- ``key:pattern:flags`` - Pattern with regex flags
- ``key:regex_pattern`` - Regular expression patterns

Supported flags:
- ``i`` - Case-insensitive matching
- ``m`` - Multiline matching
- ``s`` - Dot matches all characters including newlines

Common Examples
---------------

Filter by Gene Product
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Keep only kinase genes
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.gff3 --grep product:kinase

    # Keep transporter genes
    gfftk convert -i input.gff3 -f genome.fasta -o transporters.gff3 --grep product:transporter

    # Remove hypothetical proteins
    gfftk convert -i input.gff3 -f genome.fasta -o known_proteins.gff3 \
        --grepv product:"hypothetical.*protein"

Filter by Annotation Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Keep only augustus predictions
    gfftk convert -i input.gff3 -f genome.fasta -o augustus_only.gff3 --grep source:augustus

    # Remove augustus predictions
    gfftk convert -i input.gff3 -f genome.fasta -o no_augustus.gff3 --grepv source:augustus

    # Keep high-confidence sources
    gfftk convert -i input.gff3 -f genome.fasta -o high_conf.gff3 --grep source:genemark

Filter by Location
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Genes on specific chromosome
    gfftk convert -i input.gff3 -f genome.fasta -o chr1_genes.gff3 --grep contig:chr1

    # Genes on multiple chromosomes (regex)
    gfftk convert -i input.gff3 -f genome.fasta -o main_chrs.gff3 --grep contig:"^chr[1-5]$"

    # Plus strand genes only
    gfftk convert -i input.gff3 -f genome.fasta -o plus_strand.gff3 --grep strand:\\+

    # Minus strand genes only
    gfftk convert -i input.gff3 -f genome.fasta -o minus_strand.gff3 --grep strand:-

Advanced Filtering
------------------

Case-Insensitive Matching
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Case-insensitive search for kinases
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.gff3 --grep product:KINASE:i

    # Case-insensitive source filtering
    gfftk convert -i input.gff3 -f genome.fasta -o augustus.gff3 --grep source:AUGUSTUS:i

Regular Expression Patterns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Genes starting with specific pattern
    gfftk convert -i input.gff3 -f genome.fasta -o pattern_genes.gff3 --grep name:"^gene[0-9]+"

    # Genes with specific functional domains
    gfftk convert -i input.gff3 -f genome.fasta -o domains.gff3 \
        --grep product:"(kinase|phosphatase|transferase)"

    # Exclude ribosomal proteins
    gfftk convert -i input.gff3 -f genome.fasta -o no_ribosomal.gff3 \
        --grepv product:"ribosomal.*protein"

Multiple Filters
~~~~~~~~~~~~~~~~

You can combine multiple filters for complex selection:

.. code-block:: bash

    # Multiple grep patterns (OR logic)
    gfftk convert -i input.gff3 -f genome.fasta -o enzymes.gff3 \
        --grep product:kinase --grep product:phosphatase

    # Combined grep and grepv (AND logic)
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.gff3 \
        --grep product:kinase --grepv source:augustus

    # Complex multi-step filtering
    gfftk convert -i input.gff3 -f genome.fasta -o complex_filter.gff3 \
        --grep contig:chr1 --grep product:kinase --grepv note:pseudogene

Filter Keys
-----------

You can filter on any annotation attribute. Common keys include:

Core Attributes
~~~~~~~~~~~~~~~

- ``product`` - Gene product/function description
- ``source`` - Annotation source (augustus, genemark, etc.)
- ``name`` - Gene name or identifier
- ``contig`` - Chromosome or contig name
- ``strand`` - DNA strand (+ or -)
- ``type`` - Feature type (mRNA, CDS, etc.)

Annotation Details
~~~~~~~~~~~~~~~~~~

- ``note`` - Gene notes and comments
- ``db_xref`` - Database cross-references
- ``go_terms`` - Gene Ontology terms
- ``EC_number`` - Enzyme Commission numbers
- ``gene_synonym`` - Alternative gene names

Output Formats
--------------

Filtering works with all output formats:

.. code-block:: bash

    # Filter and convert to proteins
    gfftk convert -i input.gff3 -f genome.fasta -o kinases.faa \
        --output-format proteins --grep product:kinase

    # Filter and convert to GTF
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.gtf \
        --output-format gtf --grepv source:augustus

    # Filter and convert to TBL
    gfftk convert -i input.gff3 -f genome.fasta -o filtered.tbl \
        --output-format tbl --grep contig:chr1

    # Filter and extract transcripts
    gfftk convert -i input.gff3 -f genome.fasta -o transcripts.fasta \
        --output-format transcripts --grep product:kinase

Tips and Best Practices
-----------------------

1. **Test filters first**: Use ``--grep`` to see what matches before using ``--grepv``
2. **Quote complex patterns**: Use quotes around patterns with spaces or special characters
3. **Use anchors**: Use ``^`` and ``$`` for exact matches (e.g., ``^chr1$`` vs ``chr1``)
4. **Combine logically**: Multiple ``--grep`` = OR logic, ``--grep`` + ``--grepv`` = AND logic
5. **Validate regex**: Test complex regex patterns with online tools before using
6. **Case sensitivity**: Remember to add ``:i`` flag for case-insensitive matching

Error Handling
--------------

Common issues and solutions:

- **Invalid regex**: Check your regex syntax if you get pattern errors
- **No matches**: Verify the key name and pattern are correct
- **Case sensitivity**: Add ``:i`` flag if case doesn't match
- **Special characters**: Escape special regex characters with backslash
- **Empty results**: Check that your filter criteria aren't too restrictive

For more examples, see the :doc:`tutorial` documentation.

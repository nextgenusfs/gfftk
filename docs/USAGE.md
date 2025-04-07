# BUSCOlite Usage Guide

BUSCOlite is a streamlined implementation of BUSCO (Benchmarking Universal Single-Copy Orthologs) analysis specifically designed for genome annotation workflows. This guide explains how to use BUSCOlite for various tasks.

## Installation

### Using pip

```bash
python -m pip install buscolite
```

### From source

```bash
git clone https://github.com/nextgenusfs/buscolite.git
cd buscolite
python -m pip install .
```

## Dependencies

BUSCOlite has the following dependencies:

- [augustus](https://github.com/Gaius-Augustus/Augustus) (for genome mode)
- [miniprot](https://github.com/lh3/miniprot) (for genome mode)
- [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/index.html) (≥0.10.15)
- [pyfastx](https://github.com/lmdu/pyfastx) (≥2.0.0)
- [natsort](https://pypi.org/project/natsort/)

## Getting BUSCO Lineages

BUSCOlite requires BUSCO lineage data to run. You can download these from the BUSCO website:

- [BUSCO v5 lineages](https://busco-data.ezlab.org/v5/data/lineages/)
- [BUSCO v4 lineages](https://busco-data.ezlab.org/v4/data/lineages/)

Download the appropriate lineage for your organism, for example:

```bash
# Download fungi lineage
wget https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2020-09-10.tar.gz
# Extract
tar -xzvf fungi_odb10.2020-09-10.tar.gz
```

## Command Line Usage

### Basic Usage

```bash
buscolite -i genome.fasta -o output_name -m genome -l /path/to/fungi_odb10 -c 8
```

### Required Arguments

- `-i, --input`: Input sequence file in FASTA format (genome or proteome)
- `-o, --out`: Output name prefix for result files
- `-m, --mode`: Analysis mode, either 'genome' or 'proteins'
- `-l, --lineage`: Path to the BUSCO lineage data

### Optional Arguments

- `-c, --cpus`: Number of CPU threads to use (default: 1)
- `-s, --species`: Augustus species model to use (default: anidulans)
- `-f, --flanks`: Length of flanking regions for gene prediction (default: 2000)

### Examples

#### Genome Mode

Analyze a genome using the fungi lineage:

```bash
buscolite -i genome.fasta -o mygenome -m genome -l /path/to/fungi_odb10 -c 8 -s anidulans
```

#### Protein Mode

Analyze a proteome using the fungi lineage:

```bash
buscolite -i proteins.fasta -o myproteins -m proteins -l /path/to/fungi_odb10 -c 8
```

## Output Files

BUSCOlite generates the following output files:

- `<output_name>.buscolite.gff3`: GFF3 file with BUSCO gene annotations (genome mode only)
- `<output_name>.buscolite.tsv`: Tab-separated summary of BUSCO results
- `<output_name>.buscolite.json`: Detailed results in JSON format

### TSV Output Format

The TSV file contains the following columns:

1. BUSCO ID
2. Status (Complete, Fragmented, Missing)
3. Contig/Sequence ID
4. Start position
5. End position
6. Strand
7. Score
8. Length

### JSON Output Format

The JSON file contains detailed information about each BUSCO gene, including:

- Sequence coordinates
- Scores
- Status
- Protein sequences
- Exon/intron structure (for genome mode)

## Python API Usage

BUSCOlite can also be used as a Python library:

```python
from buscolite.busco import run_busco

results = run_busco(
    input_file="genome.fasta",
    output="output_name",
    lineage="/path/to/fungi_odb10",
    mode="genome",
    cpus=8,
    species="anidulans",
    flanks=2000
)

# Print summary
print(f"Complete: {results['complete']}")
print(f"Fragmented: {results['fragmented']}")
print(f"Missing: {results['missing']}")
print(f"Total: {results['total']}")

# Access individual BUSCO results
for busco_id, data in results['buscos'].items():
    print(f"{busco_id}: {data['status']}")
```

For more details on the Python API, see the [API documentation](API.md).

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure all dependencies are installed and in your PATH.
2. **Invalid lineage**: Make sure the lineage directory contains all required files.
3. **Augustus issues**: Some conda versions of Augustus have non-functional PPX/--proteinprofile mode.

### Checking Dependencies

```bash
# Check miniprot
miniprot --version

# Check Augustus
augustus --version

# Check Python dependencies
python -c "import pyhmmer; print(pyhmmer.__version__)"
python -c "import pyfastx; print(pyfastx.__version__)"
```

## Comparison with BUSCO

BUSCOlite is not meant to be a replacement for BUSCO. It is specifically designed for use in genome annotation pipelines like Funannotate. Key differences:

1. **Simplified workflow**: Focused on gene prediction rather than comprehensive analysis
2. **Fewer dependencies**: Requires fewer external tools
3. **Python API**: Designed to be easily integrated into other Python applications
4. **No metaeuk**: Uses miniprot and Augustus instead of metaeuk for gene prediction

For general use cases, the original [BUSCO](https://busco.ezlab.org) is recommended.

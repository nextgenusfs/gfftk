[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/gfftk.svg)](https://github.com/nextgenusfs/gfftk/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/gfftk)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/nextgenusfs/gfftk/actions/workflows/tests.yml/badge.svg)](https://github.com/nextgenusfs/gfftk/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/nextgenusfs/gfftk/branch/master/graph/badge.svg)](https://codecov.io/gh/nextgenusfs/gfftk)

# GFFtk: genome annotation tool kit

GFFtk is a comprehensive toolkit for working with genome annotation files in GFF3, GTF, and TBL formats. It provides powerful conversion, filtering, and manipulation capabilities for genomic data.

## Features

- **Format Conversion**: Convert between GFF3, GTF, TBL, and GenBank formats
- **Sequence Extraction**: Extract protein and transcript sequences from annotations
- **Advanced Filtering**: Filter annotations using flexible regex patterns
- **Consensus Models**: Generate consensus gene models from multiple sources
- **File Manipulation**: Sort, sanitize, and rename features in annotation files

## Installation

To install release versions use the pip package manager:
```bash
python -m pip install gfftk
```

To install the most updated code in master you can run:
```bash
python -m pip install git+https://github.com/nextgenusfs/gfftk.git
```

## Quick Start

### Basic Format Conversion
```bash
# Convert GFF3 to GTF
gfftk convert -i input.gff3 -f genome.fasta -o output.gtf

# Extract protein sequences
gfftk convert -i input.gff3 -f genome.fasta -o proteins.faa --output-format proteins
```

### Advanced Filtering
```bash
# Keep only kinase genes
gfftk convert -i input.gff3 -f genome.fasta -o kinases.gff3 --grep product:kinase

# Remove augustus predictions
gfftk convert -i input.gff3 -f genome.fasta -o filtered.gff3 --grepv source:augustus

# Case-insensitive filtering with regex
gfftk convert -i input.gff3 -f genome.fasta -o results.gff3 --grep product:KINASE:i

# Combined filtering
gfftk convert -i input.gff3 -f genome.fasta -o filtered.gff3 \
    --grep product:kinase --grepv source:augustus
```

### Filter Pattern Syntax
- `key:pattern` - Basic string matching
- `key:pattern:i` - Case-insensitive matching
- `key:regex` - Regular expression patterns
- Multiple `--grep` or `--grepv` flags for complex filtering

Common filter keys: `product`, `source`, `name`, `note`, `contig`, `strand`, `type`, `db_xref`, `go_terms`

For more examples and detailed documentation, see the [tutorial](docs/tutorial.rst).

## Development

### Code Formatting

This project uses [pre-commit](https://pre-commit.com/) to ensure code quality and consistency. The pre-commit hooks run Black (code formatter), isort (import sorter), and flake8 (linter).

To set up pre-commit:

1. Install pre-commit:

```bash
pip install pre-commit
```

2. Install the git hooks:

```bash
pre-commit install
```

3. (Optional) Run against all files:

```bash
pre-commit run --all-files
```

After installation, the pre-commit hooks will run automatically on each commit to ensure your code follows the project's style guidelines.

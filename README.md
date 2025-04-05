[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/gfftk.svg)](https://github.com/nextgenusfs/gfftk/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/gfftk)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/nextgenusfs/gfftk/actions/workflows/tests.yml/badge.svg)](https://github.com/nextgenusfs/gfftk/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/nextgenusfs/gfftk/branch/master/graph/badge.svg)](https://codecov.io/gh/nextgenusfs/gfftk)

# GFFtk: genome annotation tool kit


To install release versions use the pip package manager, like so:
```
python -m pip install gfftk
```

To install the most updated code in master you can run:
```
python -m pip install git+https://github.com/nextgenusfs/gfftk.git
```

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

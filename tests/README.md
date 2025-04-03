# gfftk Tests

This directory contains tests for the gfftk package.

## Running Tests

To run the tests, use pytest:

```bash
cd /Users/jon/software/gfftk
python -m pytest
```

To run tests with verbose output:

```bash
python -m pytest -v
```

To run a specific test file:

```bash
python -m pytest tests/unit/test_gff.py
```

To run a specific test:

```bash
python -m pytest tests/unit/test_gff.py::TestGFFParsing::test_detect_format
```

## Test Structure

- `unit/`: Unit tests for individual functions and classes
  - `test_basic.py`: Basic tests to verify pytest is working
  - `test_gff.py`: Tests for the gff module
  - `test_consensus.py`: Tests for the consensus module
  - `test_convert.py`: Tests for the convert module

- `data/`: Test data files
  - `sample.gff3`: Sample GFF3 file for testing
  - `sample.fasta`: Sample FASTA file for testing

## Adding Tests

When adding new tests, follow these guidelines:

1. Create a new test file in the appropriate directory (e.g., `tests/unit/test_new_module.py`)
2. Use the pytest framework for writing tests
3. Use descriptive test names that indicate what is being tested
4. Use fixtures for common test data
5. Clean up any temporary files created during tests

## Test Coverage

To run tests with coverage reporting:

```bash
python -m pytest --cov=gfftk
```

To generate a detailed HTML coverage report:

```bash
python -m pytest --cov=gfftk --cov-report=html
```

This will create a `htmlcov` directory with the coverage report.

# Contributing to GFFtk

Thank you for considering contributing to GFFtk! This document provides guidelines and instructions for contributing.

## Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```
   git clone https://github.com/your-username/gfftk.git
   cd gfftk
   ```
3. Install the package in development mode:
   ```
   pip install -e .
   ```
4. Install test dependencies:
   ```
   pip install -r test-requirements.txt
   ```

## Running Tests

We use pytest for testing. To run the tests:

```
python -m pytest
```

To run tests with coverage:

```
python -m pytest --cov=gfftk
```

To generate a coverage report:

```
python -m pytest --cov=gfftk --cov-report=html
```

This will create an HTML report in the `htmlcov` directory.

## Code Style

We follow the Black code style. To format your code:

```
black gfftk tests
```

## Pull Request Process

1. Ensure your code passes all tests
2. Update the documentation if necessary
3. Add or update tests as appropriate
4. Submit a pull request

## Continuous Integration

We use GitHub Actions for continuous integration. All pull requests will be automatically tested.

## Code of Conduct

Please be respectful and considerate of others when contributing to this project.

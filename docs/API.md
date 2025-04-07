# BUSCOlite API Documentation

BUSCOlite provides a Python API for running BUSCO analysis from within Python applications. This document describes the main functions and classes available in the API.

## Core Functions

### `buscolite.busco.run_busco`

The main function to run a BUSCO analysis.

```python
from buscolite.busco import run_busco

results = run_busco(
    input_file="genome.fasta",
    output="output_name",
    lineage="/path/to/busco_lineage",
    mode="genome",
    cpus=8,
    species="anidulans",
    flanks=2000
)
```

**Parameters:**

- `input_file` (str): Path to the input sequence file (genome or proteome)
- `output` (str): Output name prefix for result files
- `lineage` (str): Path to the BUSCO lineage data
- `mode` (str): Analysis mode, either 'genome' or 'proteins'
- `cpus` (int, optional): Number of CPU threads to use (default: 1)
- `species` (str, optional): Augustus species model to use (default: 'anidulans')
- `flanks` (int, optional): Length of flanking regions for gene prediction (default: 2000)

**Returns:**

- `dict`: Dictionary containing BUSCO results

### `buscolite.busco.check_lineage`

Verify that a BUSCO lineage directory contains all required files and directories.

```python
from buscolite.busco import check_lineage

valid, message = check_lineage("/path/to/busco_lineage")
if valid:
    print("Lineage is valid")
else:
    print(f"Lineage is invalid: {message}")
```

**Parameters:**

- `lineage` (str): Path to the BUSCO lineage directory

**Returns:**

- `tuple`: (bool, str) - Boolean indicating if the lineage is valid, and an error message if not

### `buscolite.busco.load_config`

Load the BUSCO dataset configuration file.

```python
from buscolite.busco import load_config

config = load_config("/path/to/busco_lineage")
print(f"Dataset name: {config['name']}")
```

**Parameters:**

- `lineage` (str): Path to the BUSCO lineage directory containing the dataset.cfg file

**Returns:**

- `dict`: Dictionary containing the configuration parameters from dataset.cfg

### `buscolite.busco.load_cutoffs`

Load the BUSCO score and length cutoffs from the lineage directory.

```python
from buscolite.busco import load_cutoffs

cutoffs = load_cutoffs("/path/to/busco_lineage")
for busco_id, values in cutoffs.items():
    print(f"{busco_id}: score={values['score']}, length={values['length']}")
```

**Parameters:**

- `lineage` (str): Path to the BUSCO lineage directory containing the cutoff files

**Returns:**

- `dict`: Dictionary containing the score and length cutoffs for each BUSCO model

## Search Functions

### `buscolite.search.hmmer_search`

Search multiple protein sequences against an HMM profile using pyhmmer.

```python
from buscolite.search import hmmer_search
import pyhmmer
from pyhmmer.easel import Alphabet

# Prepare sequences
alphabet = Alphabet.amino()
sequences = [seq.digitize(alphabet) for seq in protein_sequences]

results = hmmer_search("/path/to/hmm_file.hmm", sequences)
```

**Parameters:**

- `hmmfile` (str): Path to the HMM profile file
- `sequences` (list): List of digitized protein sequences to search

**Returns:**

- `list`: List of dictionaries containing search results

### `buscolite.search.hmmer_search_single`

Search a single protein sequence against an HMM profile using pyhmmer.

```python
from buscolite.search import hmmer_search_single

results = hmmer_search_single("/path/to/hmm_file.hmm", "MVNLKPTSAGRTWLKTIIIGVISAIILVVVIVIILIITSRRLNR")
```

**Parameters:**

- `hmmfile` (str): Path to the HMM profile file
- `seq` (str): Protein sequence to search

**Returns:**

- `list`: List of dictionaries containing search results

### `buscolite.search.merge_overlapping_hits`

Merge overlapping or nearby hits from BLAST or miniprot searches.

```python
from buscolite.search import merge_overlapping_hits

hits = [
    {"coords": (100, 200), "score": 10},
    {"coords": (150, 250), "score": 20},
    {"coords": (300, 400), "score": 30}
]

merged_hits = merge_overlapping_hits(hits, fluff=10000)
```

**Parameters:**

- `queryList` (list): List of dictionaries containing hit information with 'coords' key
- `fluff` (int, optional): Maximum distance between hits to consider them for merging (default: 10000)

**Returns:**

- `list`: List of merged hits

## Utility Functions

### `buscolite.utilities.runprocess`

Run a subprocess command and capture output.

```python
from buscolite.utilities import runprocess

stdout, stderr = runprocess(["ls", "-la"], cwd="/path/to/directory")
```

**Parameters:**

- `cmd` (list): Command to run as a list of strings
- `cwd` (str, optional): Working directory for the command

**Returns:**

- `tuple`: (stdout, stderr) - Standard output and standard error from the command

### `buscolite.utilities.execute`

Execute a command and yield output lines as they become available.

```python
from buscolite.utilities import execute

for line in execute(["grep", "pattern", "file.txt"]):
    print(line)
```

**Parameters:**

- `cmd` (list): Command to run as a list of strings

**Yields:**

- `str`: Lines of output from the command

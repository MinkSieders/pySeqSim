# pySeqSim

## Description
**pySeqSim** is a tool for the calculation and analysis of sequence similarity networks (SSNs) and distance matrices. It expands on the original pySSN script by providing more modular, flexible analysis of sequence relationships, clustering, and color mapping.

## Main Commands

**SSM** (`SSM`): Computes a distance matrix from input sequences using either Levenshtein or Alignment methods.

**SSN-X** (`SSN-X`): Constructs a sequence similarity network (SSN) and supports dimensionality reduction and clustering.

**C-DICT** (`C-DICT`): Generates a color dictionary for SSN visualization based on FASTA headers.

**C-MAP** (`C-MAP`): Maps sequence identifiers to color codes for SSN visualization.

## Input
Commands are accessed via a command-line interface using Python’s `argparse`. Each subcommand has specific flags and options. Use `--help` for full details.

---

## Usage | SSM

`python pySeqSim.py SSM`

Generates a distance matrix from sequences or loads a previous checkpoint.

### Required Arguments:

- `--input`, `-i` (str): Input FASTA file or distance matrix `.csv` checkpoint.
- `--metric` (str): `Alignment` or `Levenshtein`.

### Optional Arguments:

- `--threads`, `-T` (int): Number of threads to use. Default = all available CPUs.

### Example:

```bash
python pySeqSim.py SSM -i example.fasta --metric Alignment --threads 4

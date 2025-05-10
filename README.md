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


---

## Usage | SSN-X

python pySeqSim.py SSN-X

Builds and visualizes a sequence similarity network using a distance matrix or input sequences.

Required Argument:
--input, -i (str): Input file (.fasta, .csv, or .npy).

Optional Arguments:
--metric (str): "Alignment" or "Levenshtein". Only for FASTA inputs.
--threads, -T (int): Number of threads (default: CPU count).
--colors (str): TXT file with color mapping or use "cluster" for DBSCAN-based clustering.
--reducer (str): "UMAP", "tSNE", or "PCA" (default: "UMAP").
--cluster_analysis (bool): Whether to perform DBSCAN clustering (default: True).
--dbs_epsilon (float): DBSCAN epsilon (default: 0.05).
--dbs_min_sample (float): DBSCAN min samples (default: 3).
--dbs_leaf_size (float): DBSCAN leaf size (default: 30).
--faa_file (str): Original sequence file if input is a checkpoint.
--default_color (str): HEX code for default color (default: "#808080").
--legend_keys (str): Optional dictionary file to label plot legend.

Example:

python pySeqSim.py SSN-X -i example.fasta --metric Alignment --reducer UMAP --colors cluster --dbs_epsilon 0.07

---

## Usage | C-DICT

python pySeqSim.py C-DICT

Creates a color dictionary from FASTA file headers.

Required Arguments:
--input, -i (str): FASTA input file.
--method, -m (str): One of "unique_names", "starts_with", or "ends_with".

Optional Arguments:
--output, -o (str): Output TXT file name (default: "color_dictionary.txt").
--black_list, -bl (str): TXT file with sequence IDs to ignore.
--white_list, -wl (str): Required for "starts_with" or "ends_with" methods.

Example:

python pySeqSim.py C-DICT -i sequences.fasta -m starts_with -wl whitelist.txt

---

## Usage | C-MAP

python pySeqSim.py C-MAP

Maps sequences to colors using a color dictionary.

Required Arguments:
--dict, -d (str): Color dictionary TXT file.
--seqs, -s (str): Input FASTA file.

Optional Arguments:
--output, -o (str): Output TXT file (default: "color_mapping.txt").
--default_color (str): HEX fallback color (default: "#808080").

Example:

python pySeqSim.py C-MAP -d color_dictionary.txt -s sequences.fasta --default_color "#BBBBBB"

---

## Notes

- Make sure input files are correctly formatted (e.g., FASTA headers are clean and consistent).
- Distance matrices are cached to avoid redundant computation.
- All visualizations are generated using matplotlib and seaborn.

# Expansion of the pySNN script for calculation of SSNs
# Original Author: Paul Jannis Zurek
# Additional edits: Mink Sieders (m.sieders@uva.nl)

# Allows coloring by cluster & exports clustered sequences
# Additional visualization options added

#!/usr/bin/env python3

## TODO Txt annotation of selected nodes


from Bio import SeqIO
from skbio.alignment import StripedSmithWaterman
import argparse
from Levenshtein import distance as LevenDist
import pandas as pd
import seaborn as sns
import umap
from sklearn.manifold import TSNE
from collections import Counter
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import MinMaxScaler
import csv
import os
import matplotlib.colors as mcolors
from collections import defaultdict
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import itertools
from tqdm import tqdm
import multiprocessing

scaler = MinMaxScaler()

def broad_colormap_generation(n):
    colormap = plt.get_cmap("gist_rainbow")
    return ["#" + "".join(f"{int(c * 255):02X}" for c in colormap(i / n)[:3]) for i in range(n)]


def write_colors(sequences, keywords_colors, output_path):
    
    def determine_color(sequence_name, keywords_colors, generated_colors):
        default_color = hex_default
        sequence_name_lower = sequence_name.lower()[0:50]
        for keyword in keywords_colors.keys():
            if keyword == '':
                return default_color
            if keyword in sequence_name_lower:
                color = keywords_colors[keyword]
                return color if color else generated_colors.get(keyword, default_color)
        return default_color
    
    # Separate keywords that have colors from those that don't
    specified_colors = {k: v for k, v in keywords_colors.items() if v is not None}
    unspecified_keywords = [k for k, v in keywords_colors.items() if v is None]

    # Generate colors for unspecified keywords
    if unspecified_keywords:
        generated_colors = {k: mcolors.rgb2hex(c) for k, c in zip(unspecified_keywords, broad_colormap_generation(len(unspecified_keywords)))}
    else:
        generated_colors = {}

    with open(output_path, 'w') as f:
        for seq in sequences:
            seq_name = seq.split('\n')[0]
            color = determine_color(seq_name, keywords_colors, generated_colors)
            f.write(color + '\n')


def process_clusters():

    print(f"Retrieving {faa_file} for original sequences.")

    with open(cluster_file_name, "r") as f:
        clusters = f.read().splitlines()

    fasta_records = [(record.description, str(record.seq)) for record in SeqIO.parse(faa_file, "fasta")]

    if len(clusters) != len(fasta_records):
        raise ValueError("The number of clusters does not match the number of FASTA records!")

    data = [(int(cluster), header, sequence) for cluster, (header, sequence) in zip(clusters, fasta_records)]
    sorted_data = sorted(data, key=lambda x: x[0])
    output_file = f"{name}-{reducer}-clusters.tsv"

    with open(output_file, "w") as out:
        for cluster, header, sequence in sorted_data:
            out.write(f"{cluster}\t{header}\t{sequence}\n")

    print(f"Intermediate cluster file with sequences written to {output_file}. Run {output_file} as input to generate "
          f".faa files for all clusters individually")


def extract_all_clusters_to_fasta():

    output_dir = f"{name}"
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        clusters = {}

        with open(input_file, "r") as tsv:
            reader = csv.reader(tsv, delimiter="\t")

            for row in reader:
                if len(row) < 3:
                    continue

                cluster, header, sequence = row
                if cluster not in clusters:
                    clusters[cluster] = []
                clusters[cluster].append((header, sequence))

        for cluster, sequences in clusters.items():
            num_cluster_for_formatting = f"{int(cluster):05d}"
            output_fasta = os.path.join(output_dir, f"cluster_{num_cluster_for_formatting}.faa")
            with open(output_fasta, "w") as fasta:
                for header, sequence in sequences:
                    fasta.write(f">{header}\n{sequence}\n")
            print(f"Extracted {len(sequences)} sequences for cluster {cluster} to {output_fasta}")

    except FileNotFoundError:
        print(f"File {input_file} not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def all_scores(query_nr, record):   # Calculates all scores for distance_matrix
    global query_lst
    aln = query_lst[query_nr](record)
    score = aln.optimal_alignment_score
    query_length = len(aln.query_sequence)
    aligned_query_length = aln.query_end - aln.query_begin + 1
    coverage = aligned_query_length / query_length
    aln_query = aln.aligned_query_sequence
    aln_target = aln.aligned_target_sequence
    aln_length = len(aln_query)
    same_aa = sum(e1 == e2 for e1, e2 in zip(aln_query, aln_target))
    ident = same_aa / aln_length
    return [score, ident, coverage]


def cluster_and_save(embedding):
    clustering = DBSCAN(eps=0.05, min_samples=3,leaf_size=30).fit(embedding)
    cluster_labels = clustering.labels_
    with open(cluster_file_name, "w") as f:
        for label in cluster_labels:
            f.write(f"{label}\n")

    return cluster_labels


def plot_scatter_colored(embedding):
    
    def darken_color(hex_color, factor=0.8):

        r = int(hex_color[1:3], 16)
        g = int(hex_color[3:5], 16)
        b = int(hex_color[5:7], 16)

        r = int(r * factor)
        g = int(g * factor)
        b = int(b * factor)

        return '#{:02x}{:02x}{:02x}'.format(r, g, b)
    
    if perform_clustering:
        clusters = cluster_and_save(embedding)
    else:
        clusters = None

    if legend_keys is not None:
        legend_file = legend_keys
        legend = {}
        with open(legend_file, 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    key, value = line.split(',')
                    legend[key] = value

    if grouping == "cluster":
        if clusters is None:
            raise ValueError("Grouping set to cluster while --cluster_analysis is put to False.")

        unique_clusters = sorted(set(clusters))
        cluster_color_map = {}
        cmap = plt.get_cmap('tab20')
        num_colors = len(unique_clusters) - (1 if -1 in unique_clusters else 0)
        color_iter = iter(cmap(np.linspace(0, 1, num_colors)))
        for cluster in unique_clusters:
            if cluster == -1:
                cluster_color_map[cluster] = hex_default
            else:
                rgba = next(color_iter)
                hex_color = '#{:02x}{:02x}{:02x}'.format(
                    int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
                )
                cluster_color_map[cluster] = hex_color

        colors = [cluster_color_map[c] for c in clusters]

    else:
        colors = ["k" for _ in range(len(embedding[:, 1]))]
        if grouping is not None:
            with open(grouping, 'r') as f:
                lines = f.readlines()
            if len(lines) != len(colors):
                raise ValueError("Number of groupings does not match number of embedded sequences.")
            colors = [l.strip("\n") for l in lines]
            try:
                colors = [int(c) for c in colors]
            except Exception:
                pass

    color_counts = Counter(colors)
    sorted_colors = sorted(color_counts, key=color_counts.get)
    most_frequent_color = sorted_colors[-1]

    if grouping == "cluster":
        plt.figure()
        for cluster in unique_clusters:
            mask = np.array(clusters) == cluster
            if cluster == -1:
                plt.scatter(
                    embedding[mask, 0], embedding[mask, 1], c=cluster_color_map[cluster],
                    label=f"Noise (Cluster {cluster})", s=1, alpha=0.4
                )
            else:
                plt.scatter(
                    embedding[mask, 0], embedding[mask, 1], c=cluster_color_map[cluster],
                    label=f"Cluster {cluster}", s=4, alpha=0.75
                )

                centroid = embedding[mask].mean(axis=0)

                plt.text(
                    centroid[0], centroid[1], str(cluster),
                    color=darken_color(cluster_color_map[cluster], factor=0.8),
                    fontsize=10, fontweight="bold", ha="center", va="center"
                )
    else:
        plt.figure()
        for color_i in reversed(sorted_colors):
            mask = np.array(colors) == color_i
            print(f"Color: {color_i}, Mask sum: {np.sum(mask)}")
            alpha = 1 if color_i != most_frequent_color else 0.2
            size = 4 if color_i != most_frequent_color else 1
            plt.scatter(embedding[mask, 0], embedding[mask, 1], c=[color_i], s=size, alpha=alpha, edgecolors=None)

    plt.xticks([])
    plt.yticks([])

    if grouping != "cluster":
        if legend_keys is not None:
            handles = []
            for key, value in legend.items():
                handles.append(plt.scatter([], [], color=value, label=key, marker='o'))
            plt.legend(handles=handles, loc="upper left", bbox_to_anchor=(1, 1), ncol=(len(handles) // 20 + 1))

    plt.savefig(f"{name}-{reducer}.svg", bbox_inches="tight", dpi=600)


def distance_matrix(records):
    rec_lst = records.copy()
    global query_lst
    query_lst = [StripedSmithWaterman(rec) for rec in rec_lst]
    pool = multiprocessing.Pool(processes=threads)
    score_lst_lst = []
    N_rec = len(rec_lst)
    for i in tqdm(range(N_rec)):
        rec_lst.pop(0)
        score_lst = pool.starmap(all_scores, zip(itertools.repeat(i), rec_lst))
        identlst = [1-elem[1] for elem in score_lst] # 1-identity = dissimilarity
        score_lst_lst.append(identlst)
    pool.close()
    print('finished generating the alignment score matrix')
    return score_lst_lst


def convert_distance_matrix(sparse_matrix):   
    elem = len(sparse_matrix)
    full_matrix = [[0 for i in range(elem)] for j in range(elem)]
    for i in range(len(sparse_matrix)):
        for j in range(len(sparse_matrix[i])):
            full_matrix[i][j+1+i] = sparse_matrix[i][j]
            full_matrix[j+1+i][i] = sparse_matrix[i][j]
    return full_matrix


def levenhstein_distance_matrix(records):
    rec_lst = records.copy()
    score_lst_lst = []
    N_rec = len(rec_lst)
    pool = multiprocessing.Pool(processes=threads)
    max_dist = 0
    for i in tqdm(range(N_rec)):
        ref = rec_lst.pop(0)
        score_lst = pool.starmap(LevenDist, zip(itertools.repeat(ref), rec_lst))
        m = max(score_lst, default=0)
        if m > max_dist:
            max_dist = m
        score_lst_lst.append(score_lst)
    pool.close()
    # Scale 0-1
    for i in range(len(score_lst_lst)):
        for j in range(len(score_lst_lst[i])):
            score_lst_lst[i][j] = score_lst_lst[i][j] / max_dist
    print('finished generating the levenshtein score matrix')
    return score_lst_lst


def calc_reduction(DM):
    if reducer == "UMAP":
        model = umap.UMAP(metric="precomputed")
    elif reducer == "tSNE":
        model = TSNE(n_components=2)
    elif reducer == "PCA":
        model = PCA(n_components=2)

    embedding = model.fit_transform(DM)
    scaler.fit(embedding)
    norm_embedding = scaler.transform(embedding)
    np.save(f"{name}-{reducer}checkpoint.npy", norm_embedding, allow_pickle=True)

    return norm_embedding


def color_dict_generation(in_file, out_file, black_list_file, method, white_list_file=None):
    unique_names = set()

    with open(black_list_file, 'r') as f:
        black_list = {line.strip().lower() for line in f}

    if method in ('starts_with', 'ends_with'):
        if not white_list_file:
            raise ValueError(f"Whitelist file must be provided when using method '{method}'.")
        with open(white_list_file, 'r') as f:
            white_list = {line.strip().lower() for line in f}
    else:
        white_list = set()

    with open(in_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                name_lc = name.lower()

                if any(bl in name_lc for bl in black_list):
                    continue

                if method == 'starts_with':
                    if not any(name_lc.startswith(wl) for wl in white_list):
                        continue
                elif method == 'ends_with':
                    if not any(name_lc.endswith(wl) for wl in white_list):
                        continue
                elif method != 'unique_names':
                    raise ValueError("Invalid method. Choose from ['unique_names', 'starts_with', 'ends_with']")

                unique_names.add(name)

    color_list = broad_colormap_generation(len(unique_names))
    assigned_colors = dict(zip(unique_names, color_list))

    with open(out_file, 'w') as f:
        for name, color in assigned_colors.items():
            f.write(f"{name},{color}\n")


def sequence_similarity_matrix(DM, names):

    # Revert dissimilarity to similarity matrix

    sim_DM = 1 - DM

    plt.figure(figsize=(10, 8))
    plt.imshow(sim_DM, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Sequence Identity')
    plt.xticks([])
    plt.yticks([])
    plt.xticks(np.arange(len(names)), names, rotation=90, fontsize=7)
    plt.yticks(np.arange(len(names)), names, fontsize=7)
    plt.tight_layout()
    file_name_ssm = f"{name}-{metric}-SSM.svg"
    plt.savefig(file_name_ssm)

    print(f"Similarity matrix saved as '{file_name_ssm}'")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""pySeqSim [Calculation & Analysis of SSNs and Distance Matrices].
                                     Expansion of pySSN script.
                                     Original pySSN Author: Paul Zurek.
                                     Author pySeqSim: Mink Sieders (m.sieders@uva.nl).
                                     Version 1.0""")

    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_ssm = subparsers.add_parser("SSM", help= "Run the sequence similarity matrix analysis")
    parser_ssm.add_argument('-T', '--threads', type=int, default=0,
                            help='Number of threads to execute in parallel. Defaults to CPU count.')
    parser_ssm.add_argument('-i', '--input', type=str, help="""Please provide one of the following input files:\n 
                                                              FASTA: List of records to calculate distance matrix from.\n
                                                              CSV: Distance matrix checkpoint.
                                                              """, required=True)
    parser_ssm.add_argument('--metric', type=str, choices=['Levenshtein','Alignment'], default='Alignment', help='Metic used for distance calculation: Levenshtein or Alignment. Use Levenshtein for close sequences and Alignment for less homologous sequences.')

    parser_ssn = subparsers.add_parser("SSN-X", help="Run the sequence similarity network analysis")

    parser_ssn.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
    parser_ssn.add_argument('-i', '--input', type=str, help="""Please provide one of the following input files:\n 
                                                          FASTA: List of records to calculate distance matrix from.\n
                                                          CSV: Distance matrix checkpoint.
                                                          NPY: Reducer embeddings checkpoint.
                                                          """, required=True)
    parser_ssn.add_argument('--metric', type=str, choices=['Levenshtein','Alignment'], default='Alignment', help='Metic used for distance calculation: Levenshtein or Alignment. Use Levenshtein for close sequences and Alignment for less homologous sequences.')
    parser_ssn.add_argument('--colors', type=str, help='TXT file for color information. Can be used to color the SSN. If option "cluster" is given, color will be determined by DBSCAN clustering analysis.')
    parser_ssn.add_argument('--reducer', type=str, choices=['UMAP','tSNE', 'PCA'], default="UMAP", help='Choice of dimensionality reduction method: UMAP, PCA, or tSNE. Defaults to UMAP.')
    parser_ssn.add_argument('--cluster_analysis', type=bool, default=True, help='Perform clustering on SSN using DBSCAN. Does not automatically color plot according to clusters, use "--grouping cluster" for this. ')
    parser_ssn.add_argument('--dbs_epsilon', type=float, default=0.05, help='DBSCAN Epsilon setting, only relevant for clustering analysis.')
    parser_ssn.add_argument('--dbs_min_sample', type=float, default=3, help='DBSCAN Minimum Samples setting, only relevant for clustering analysis.')
    parser_ssn.add_argument('--dbs_leaf_size', type=float, default=30, help='DBSCAN Leaf Size setting, only relevant for clustering analysis.')
    parser_ssn.add_argument('--faa_file', type=str, default=None,help='Specify .fasta file where original sequences can be found on which clustering was performed. Automatically defaults to "{original_name}.faa".')
    parser_ssn.add_argument('--default_color', type=str,
                                 help="HEX Code for the default bead color in the SSN plot [Default = Grey].",
                                 default='#808080',
                                 required=False)
    parser_ssn.add_argument('--legend_keys', type=str,
                            help="Colormap dictionary file for legend keys to generate legend in the plot.",
                            default=None,
                            required=False)

    parser_keygen = subparsers.add_parser("C-DICT", help="Color dictionary generator for SNN-X plot color mapping")
    parser_keygen.add_argument('-i', '--input', type=str, default=None, help=".fasta input file for C-DICT.", required=True)
    parser_keygen.add_argument('-m', '--method', type=str, choices=['unique_names', 'starts_with', 'ends_with'],help="Method of writing a color dictionary based"
                                                                "on fasta sequence IDs", required=True)
    parser_keygen.add_argument('-o', '--output', type=str,
                               help="Output name (.txt) for color dictionary file.", default = 'color_dictionary.txt', required=False)
    parser_keygen.add_argument('-bl', '--black_list', type=str,
                               help="File path (.txt) for blacklisted sequence ID's (omitted in C-DICT)",
                               required=False)
    parser_keygen.add_argument('-wl', '--white_list', type=str,
                               help="File path (.txt) for whitelisted sequence ID's, required when using method starts_with and eds_width in the --method flag.",
                               required=False)

    parser_colormap = subparsers.add_parser("C-MAP", help="Map sequences to colors for SSN-X plots")
    parser_colormap.add_argument('-d', '--dict', type=str, help="Input dictionary file (.txt) for ColorMap.", required=True)
    parser_colormap.add_argument('-s', '--seqs', type=str, help="Input sequences file (.fasta) for ColorMap.", required=True)
    parser_colormap.add_argument('-o', '--output', type=str,
                               help="Output name (.txt) for color mapping file.", default='color_mapping.txt',
                               required=False)
    parser_colormap.add_argument('--default_color', type=str,
                                 help="HEX Code for the default bead color in the SSN-X plot [Default = Grey].", default='#808080',
                                 required=False)


    args = parser.parse_args()


    if args.command == "SSN-X":

        threads = args.threads
        metric = args.metric
        grouping = args.colors
        reducer = args.reducer
        perform_clustering = args.cluster_analysis
        dbs_epsilon = args.dbs_epsilon
        dbs_min_sample = args.dbs_min_sample
        dbs_leaf_size = args.dbs_leaf_size
        faa_file = args.faa_file
        input_file = args.input
        hex_default = args.default_color
        legend_keys = args.legend_keys

        name = ".".join(input_file.split(".")[:-1])
        ftype = input_file.split(".")[-1].lower()
        cluster_file_name = f"{name}-{reducer}-clusters.txt"
        sns.set(style='white', context='notebook', rc={'figure.figsize': (6, 6)})
        if threads == 0:
            threads = multiprocessing.cpu_count()

        print(f"Running SNN-X with input {input_file} and reducer {reducer} and {threads} threads.")

        if ftype == "fasta":
            if metric is None:
                raise NameError("Please specify a distance metric.")
            records = [str(rec.seq) for rec in SeqIO.parse(input_file, "fasta")]
            print(f"{len(records)} records loaded.")
            print(f"Calculating distance matrix via {metric}:")
            if metric == 'Alignment':
                DM = distance_matrix(records)
            elif metric == 'Levenshtein':
                DM = levenhstein_distance_matrix(records)
            fullDM = convert_distance_matrix(DM)
            fullDM = pd.DataFrame(fullDM)

            dm_file = f"{name}-{metric}DM-checkpoint.csv"
            fullDM.to_csv(dm_file)
            sequences = list(SeqIO.parse(input_file, "fasta"))
            sequence_names = [seq.id for seq in sequences]
            seqname_file = f"{dm_file}.seqnames"
            with open(seqname_file, "w") as f:
                for burlp in sequence_names:
                    f.write(f"{burlp}\n")

            embeddings = calc_reduction(fullDM)
            cluster_file_name = f"{name}-{reducer}-clusters.txt"
            plot_scatter_colored(embeddings)
            faa_file = input_file
            process_clusters()

        elif ftype == "csv":  # From DM checkpoint
            print("Loading precomputed distance matrix")
            fullDM = pd.read_csv(input_file, index_col=0)
            print("calculating embeddings")
            embeddings = calc_reduction(fullDM)
            plot_scatter_colored(embeddings)

            if faa_file is None:
                faa_file = f"{name}.fasta"
                print(f"No --faa_file provided for original sequences, assuming sequences are found in {faa_file}.")
            process_clusters()

        elif ftype == "npy":  # From embeddings checkpoint (e.g. to re-color without re-calculation OR re-perform cluster analysis)
            print("Loading precomputed embedding")
            embeddings = np.load(input_file, allow_pickle=True)
            plot_scatter_colored(embeddings)

            if faa_file is None:
                faa_file = f"{name}.fasta"
                print(f"No --faa_file provided for original sequences, assuming sequences are found in {faa_file}.")
            process_clusters()

        elif ftype == "tsv": #Processed cluster intermediate file
            print("Loading precomputed clusters determined in previous analysis.")
            extract_all_clusters_to_fasta()

        else:
            print("Unrecognized input file.")


    elif args.command == "C-DICT":
        out_file = args.output
        in_file = args.input
        method = args.method
        black_list = args.black_list
        white_list = args.white_list

        color_dict_generation(in_file, out_file, black_list, method, white_list)

        print(f"Color dictionary written to {out_file}")


    elif args.command == "C-MAP":

        out_map = args.output
        in_seqs = args.seqs
        in_dict = args.dict
        hex_default = args.default_color

        sequences = []
        with open(in_seqs, 'r') as f:
            current_seq = ''
            for line in f:
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(current_seq)
                    current_seq = line.strip()
                else:
                    current_seq += line.strip()
            if current_seq:
                sequences.append(current_seq)

        keywords_colors = defaultdict(lambda: None)
        with open(in_dict, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                keyword = parts[0].lower()
                color = parts[1] if len(parts) > 1 else None
                keywords_colors[keyword] = color

        write_colors(sequences, keywords_colors, out_map)

        print(f"Color mappings written to {out_map}")


    elif args.command == "SSM":
        threads = args.threads
        input_file = args.input
        metric = args.metric

        if threads == 0:
            threads = multiprocessing.cpu_count()

        name = ".".join(input_file.split(".")[:-1])
        ftype = input_file.split(".")[-1].lower()

        if ftype == "fasta":
            if metric is None:
                raise NameError("Please specify a distance metric.")
            records = [str(rec.seq) for rec in SeqIO.parse(input_file, "fasta")]
            print(f"{len(records)} records loaded.")
            print(f"Calculating distance matrix via {metric}:")
            if metric == 'Alignment':
                DM = distance_matrix(records)
            elif metric == 'Levenshtein':
                DM = levenhstein_distance_matrix(records)
            fullDM = convert_distance_matrix(DM)
            fullDM = pd.DataFrame(fullDM)

            # Save DM
            dm_file = f"{name}-{metric}DM-checkpoint.csv"
            fullDM.to_csv(dm_file)
            sequences = list(SeqIO.parse(input_file, "fasta"))
            sequence_names = [seq.id for seq in sequences]
            seqname_file = f"{dm_file}.seqnames"
            with open(seqname_file, "w") as f:
                for burlp in sequence_names:
                    f.write(f"{burlp}\n")

            # Do SSM Stuff
            sequence_similarity_matrix(fullDM, sequence_names)


        elif ftype == "csv":  # From DM checkpoint
            print("Loading precomputed distance matrix")
            fullDM = pd.read_csv(input_file, index_col=0)

            # Try to load sequence names from corresponding file
            seqname_file = f"{input_file}.seqnames"
            try:
                with open(seqname_file, "r") as f:
                    sequence_names = [line.strip() for line in f]
            except FileNotFoundError:
                raise FileNotFoundError(f"Sequence names file '{seqname_file}' not found. Make sure it was generated.")

            # Do SSM Stuff
            sequence_similarity_matrix(fullDM, sequence_names)

        else:
            print("Unrecognized input file.")

import argparse
from Bio import SeqIO
from itertools import combinations
from joblib import Parallel, delayed
import json
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import statistics
import sys
from tqdm import tqdm

def read_metadata(input_json):
    with open(input_json, 'r') as f:
        metadata = json.load(f)
    return metadata

def rolling_hash(seq, k, base=4):
    hash_value = 0
    for i in range(k):
        hash_value = hash_value * base + ord(seq[i])
    return hash_value

def update_rolling_hash(prev_hash, prev_kmer, next_kmer, k, base=4):
    return (prev_hash - ord(prev_kmer[0]) * (base ** (k - 1))) * base + ord(next_kmer[-1])

def calculate_minimizers(seqRecord, k, w):
    rev_comp_sequence = str(seqRecord.reverse_complement().seq)
    sequence = str(seqRecord.seq)

    # Initialize hashes for the first k-mer in both strands
    hash_forward = rolling_hash(sequence[:k], k)
    hash_rev_comp = rolling_hash(rev_comp_sequence[:k], k)

    hashes_forward = [hash_forward]
    hashes_rev_comp = [hash_rev_comp]

    # Calculate rolling hashes for the rest of the k-mers in both strands
    for i in range(1, len(sequence) - k + 1):
        hash_forward = update_rolling_hash(hash_forward, sequence[i-1:i-1+k], sequence[i:i+k], k)
        hash_rev_comp = update_rolling_hash(hash_rev_comp, rev_comp_sequence[i-1:i-1+k], rev_comp_sequence[i:i+k], k)
        hashes_forward.append(hash_forward)
        hashes_rev_comp.append(hash_rev_comp)

    # Find minimizers in windows of w consecutive k-mers
    minimizers = set([min(min(hashes_forward[i:i+w]), min(hashes_rev_comp[i:i+w])) for i in range(len(hashes_forward) - w + 1)])
    assert not len(minimizers) == 0, "Your chosen windowSize is larger than your shortest sequence. Reduce the windowSize and rerun the tool."
    return (seqRecord.id, minimizers)

def calculate_jaccard_distance(chunk,
                            shorter):
    distances = []
    for c in tqdm(chunk):
        seq1 = c[0]
        seq2 = c[1]
        len_intersection = len(seq1[1].intersection(seq2[1]))
        if not shorter:
            dist = float(len_intersection / (len(seq1[1]) + len(seq2[1]) - len_intersection))
        else:
            shorter = seq1[1] if len(seq1[1]) < len(seq2[1]) else seq2[1]
            dist = float(len_intersection / len(shorter))
        distances.append((seq1[0], seq2[0], dist))
    return distances

def write_distance_gml(sequence_records,
                    mash_df,
                    metadata,
                    output_gml):
    G = nx.Graph()
    sys.stderr.write("\t\tWriting nodes...\n")
    for seq in tqdm(sequence_records):
        G.add_node(seq.id, cluster=metadata[seq.id])
    sys.stderr.write("\t\tWriting edges...\n")
    for idx, row in tqdm(mash_df.iterrows()):
        G.add_edge(row['seq1'], row['seq2'], weight=row['identity'])
    nx.write_gml(G, output_gml)

def calculate_cluster_stats(cluster_df, metadata):
    seq1 = [metadata[s] for s in list(cluster_df["seq1"])]
    seq2 = [metadata[s] for s in list(cluster_df["seq2"])]
    cluster_df["seq1_cluster"] = seq1
    cluster_df["seq2_cluster"] = seq2
    cluster_df["Same cluster"] = cluster_df["seq1"] == cluster_df["seq2"]

    distances = list(cluster_df["identity"])
    clustered = [str(seq1[i]) == str(seq2[i]) for i in range(len(seq1))]
    cluster1 = list(cluster_df["seq1_cluster"])
    cluster2 = list(cluster_df["seq2_cluster"])

    withinClusterIdentities = {}
    betweenClusterIdentities = {}

    for i in tqdm(range(len(clustered))):
        if not cluster1[i] in betweenClusterIdentities:
            betweenClusterIdentities[cluster1[i]] = []
        if not cluster1[i] in withinClusterIdentities:
            withinClusterIdentities[cluster1[i]] = []
        if not cluster2[i] in betweenClusterIdentities:
            betweenClusterIdentities[cluster2[i]] = []
        if not cluster2[i] in withinClusterIdentities:
            withinClusterIdentities[cluster2[i]] = []
        if clustered[i]:
            withinClusterIdentities[cluster1[i]].append(distances[i])
            withinClusterIdentities[cluster2[i]].append(distances[i])
        else:
            betweenClusterIdentities[cluster1[i]].append(distances[i])
            betweenClusterIdentities[cluster2[i]].append(distances[i])

    means = {"within_cluster": [], "between_cluster": []}
    median = {"within_cluster": [], "between_cluster": []}
    mode = {"within_cluster": [], "between_cluster": []}
    ranges = {"within_cluster": [], "between_cluster": []}
    for clus in tqdm(withinClusterIdentities):
        if not withinClusterIdentities[clus] == []:
            means["within_cluster"].append(statistics.mean(withinClusterIdentities[clus]))
            median["within_cluster"].append(statistics.median(withinClusterIdentities[clus]))
            mode["within_cluster"].append(statistics.mode(withinClusterIdentities[clus]))
            ranges["within_cluster"].append(max(withinClusterIdentities[clus]) - min(withinClusterIdentities[clus]))
        else:
            means["within_cluster"].append(1)
            median["within_cluster"].append(1)
            mode["within_cluster"].append(1)
            ranges["within_cluster"].append(0)
    for clus in tqdm(betweenClusterIdentities):
        if not betweenClusterIdentities[clus] == []:
            means["between_cluster"].append(statistics.mean(betweenClusterIdentities[clus]))
            median["between_cluster"].append(statistics.median(betweenClusterIdentities[clus]))
            mode["between_cluster"].append(statistics.mode(betweenClusterIdentities[clus]))
            ranges["between_cluster"].append(max(betweenClusterIdentities[clus]) - min(betweenClusterIdentities[clus]))
        else:
            means["between_cluster"].append(0)
            median["between_cluster"].append(0)
            mode["between_cluster"].append(0)
            ranges["between_cluster"].append(0)
    return means, median, mode, ranges

def create_jointplot(means, median, mode, ranges,
                    clusterPlot):
    # Set the aesthetics for the plot
    sns.set(style="white", color_codes=True)
    for i in [("mean", means), ("median", median), ("mode", mode), ("range", ranges)]:
        # Create the jointplot
        jointplot = sns.jointplot(x=i[1]["within_cluster"],
                                y=i[1]["between_cluster"],
                                kind='kde')
        jointplot.ax_marg_x.set_xlim(0, 1)
        jointplot.ax_marg_y.set_ylim(0, 1)
        jointplot.ax_joint.set_xlabel('Within cluster sequence identity')
        jointplot.ax_joint.set_ylabel('Between cluster sequence identity')
        jointplot.figure.tight_layout()
        jointplot.figure.savefig(clusterPlot + "." + i[0] + ".png",
                                dpi=300)

def get_options():
    parser = argparse.ArgumentParser(description='Create visualisations of approximate between and within cluster nucleotide identities for short sequences.')
    parser.add_argument('input_fasta', help='Input FASTA file of all sequences.')
    parser.add_argument('input_json', help='Input JSON file with cluster assignments ({<sequence header>: <cluster assignment>}).')
    parser.add_argument('--clusterGML', default=None, help='Output path of GML clustering file to view with Cytoscape or similar.')
    parser.add_argument('--distanceTable', default=None, help='Output path of CSV of identities (may take a long time).')
    parser.add_argument('--clusterPlot', default=None, help='Output path of jointplot to visualise between and within cluster identities.')
    parser.add_argument('--kmerSize', type=int, default=9, help='Kmer size (default: 9).')
    parser.add_argument('--windowSize', type=int, default=20, help='Minimiser window size (default: 20).')
    parser.add_argument('--threshold', type=float, default=0.9, help='Jaccard similarity threshold (default: 0.9).')
    parser.add_argument('--threads', type=int, default=1, help='Threads for sketching and jaccard distance calculations (default: 1).')
    parser.add_argument('--shorter', dest='shorter', action='store_true', default=False,
                        help='Assess identity relative to the shorter sequence.')
    args = parser.parse_args()
    return args

def main():
    # parse command line arguements
    args = get_options()
    # read cluster assignment metadata
    metadata = read_metadata(args.input_json)
    # read fasta files
    sys.stderr.write("Reading sequences\n")
    sequence_records = list(SeqIO.parse(args.input_fasta, "fasta"))
    # extract minimisers from sequences
    sys.stderr.write(f"Extracting minimisers using {str(args.threads)} threads...\n")
    kmer_sets = Parallel(n_jobs=args.threads)(delayed(calculate_minimizers)(seq,
                                                                        args.kmerSize,
                                                                        args.windowSize) for seq in tqdm(sequence_records))
    # Calculate distances
    sys.stderr.write(f"Calculating approximate jaccard distances using {str(args.threads)} threads...\n")
    distances = []
    # Get all unique pairs of sequences
    all_comparisons = [c for c in combinations(kmer_sets, 2)]
    # Split pairs into chunks based on the number of threads
    chunks = np.array_split(all_comparisons, args.threads)
    # Calculate approximate jaccard distances
    chunk_distances = Parallel(n_jobs=args.threads)(delayed(calculate_jaccard_distance)(c, args.shorter) for c in chunks)
    for thread_distances in chunk_distances:
        distances += thread_distances
    # convert pairwise distances to dataframe
    cluster_df = pd.DataFrame(distances, columns=['seq1', 'seq2', 'identity'])
    if args.clusterPlot:
        # separate between cluster and within cluster distances
        sys.stderr.write("Making jointplot of all pairwise identities...\n")
        means, median, mode, ranges = calculate_cluster_stats(cluster_df, metadata)
        # plot jointplots of between and within cluster distances
        create_jointplot(means, median, mode, ranges,
                        args.clusterPlot)
    # filter out pairs with distances below the threshold
    if args.distanceTable or args.clusterGML:
        sys.stderr.write("Filtering out pairs of sequences between the threshold...\n")
        distances = [dist for dist in tqdm(distances) if dist[2] >= args.threshold]
    if args.distanceTable:
        # Write filtered CSV output
        sys.stderr.write("Writing distance CSV\n")
        cluster_df = pd.DataFrame(distances, columns=['seq1', 'seq2', 'identity'])
        cluster_df.to_csv(args.distanceTable, index=False)
    # Write GML output
    if args.clusterGML:
        sys.stderr.write("Writing distance GML\n")
        write_distance_gml(sequence_records,
                        cluster_df,
                        metadata,
                        args.clusterGML)

if __name__ == "__main__":
    main()
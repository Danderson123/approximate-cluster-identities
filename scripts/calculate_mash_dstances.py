import os
from Bio import SeqIO
import argparse
from approximate_cluster_identities.__main__ import read_metadata, write_distance_gml, create_jointplot, calculate_cluster_stats
from joblib import Parallel, delayed
from tempfile import TemporaryDirectory
from tqdm import tqdm
import sys
import pandas as pd
import subprocess
from itertools import combinations


def write_sequence_to_temp_file(sequence, temp_dir):
    sequenceFileName = os.path.join(temp_dir, sequence.id + ".fasta")
    SeqIO.write(sequence, sequenceFileName, "fasta")
    return sequenceFileName

def run_mash_sketch(sequence_files_path, threads, mash_path, kmer_size):
    sketch_file = os.path.join(os.path.dirname("./."), "reference.msh")
    os.system(f"{mash_path} sketch -k {kmer_size} -p {threads} -l {sequence_files_path} -o {sketch_file}")
    return sketch_file

def run_mash_dist(sketch_file, mash_path, threads, kmerSize):
    mash_output = subprocess.check_output([f"{mash_path}", "dist", "-s", "100", "-p", f"{threads}", f"{sketch_file}", f"{sketch_file}"]).decode()
    # Parse the output into a dataframe
    data = []
    lines = mash_output.split('\n')
    for line in lines:
        if line:
            elements = line.split('\t')
            data.append({
                "seq1": os.path.basename(elements[0]).replace(".fasta", ""),
                "seq2": os.path.basename(elements[1]).replace(".fasta", ""),
                "identity": 1 - float(elements[2]),  # Convert distance to identity
            })
    cluster_df = pd.DataFrame(data, columns=['seq1', 'seq2', 'identity'])
    return cluster_df

def get_options():
    parser = argparse.ArgumentParser(description='Create GML file from FASTA sequences using Mash distances.')
    parser.add_argument('input_fasta', help='Input FASTA file.')
    parser.add_argument('input_json', help='Input JSON file with metadata.')
    parser.add_argument('mash_path', help='Path to mash binary.')
    parser.add_argument('--clusterGML', default=None, help='Output path of GML clustering file to view with Cytoscape or similar.')
    parser.add_argument('--distanceTable', default=None, help='Output path of CSV of identities (may take a long time).')
    parser.add_argument('--clusterPlot', default=None, help='Output path of jointplot to visualise between and within cluster identities.')
    parser.add_argument('--kmerSize', type=int, default=11, help='Kmer size (default: 11).')
    parser.add_argument('--windowSize', type=int, default=100, help='Minimiser window size (default: 100).')
    parser.add_argument('--threshold', type=float, default=0.9, help='Jaccard similarity threshold (default: 0.9).')
    parser.add_argument('--threads', type=int, default=1, help='Threads for sketching and jaccard distance calculations (default: 1).')
    args = parser.parse_args()
    return args

def main():
    args = get_options()
    #read cluster assignment metadata
    metadata = read_metadata(args.input_json)
    # read fasta files
    sys.stderr.write("Reading sequences\n")
    sequence_records = list(SeqIO.parse(args.input_fasta, "fasta"))
    # Sketch sequences using Mash
    with TemporaryDirectory() as temp_dir:
        sequence_files = Parallel(n_jobs=args.threads)(delayed(write_sequence_to_temp_file)(seq, temp_dir) for seq in tqdm(sequence_records))
        # Write sequence file paths to a text file
        sequence_files_path = os.path.join(temp_dir, "sequence_files.txt")
        with open(sequence_files_path, 'w') as f:
            for seq_file in sequence_files:
                f.write(seq_file + '\n')
        # Sketch sequences using Mash
        sys.stderr.write("Running mash dist\n")
        # Calculate distances
        sketch_file = run_mash_sketch(sequence_files_path, args.threads, args.mash_path, args.kmerSize)
        cluster_df = run_mash_dist(sketch_file, args.mash_path, args.threads, args.kmerSize)
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
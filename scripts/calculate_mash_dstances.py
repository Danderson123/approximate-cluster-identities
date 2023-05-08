import os
from Bio import SeqIO
import argparse
from approximate_cluster_identities.__main__ import read_metadata, write_distance_gml
from joblib import Parallel, delayed
from tempfile import TemporaryDirectory
from tqdm import tqdm
import sys
import pandas as pd
from itertools import combinations


def write_sequence_to_temp_file(sequence, temp_dir):
    sequenceFileName = os.path.join(temp_dir, sequence.id + ".fasta")
    SeqIO.write(sequence, sequenceFileName, "fasta")
    return sequenceFileName

def run_mash_sketch(sequence_files_path, threads, mash_path):
    sketch_file = os.path.join(os.path.dirname("./."), "reference.msh")
    os.system(f"{mash_path} sketch -k 11 -r -p {threads} -l {sequence_files_path} -o {sketch_file}")
    return sketch_file

def run_mash_dist(file, sketch_file, mash_path):
    cmd = f"{mash_path} dist {sketch_file} {file}"
    result = os.popen(cmd).read().split("\t")
    return (os.path.basename(result[0]).replace(".fasta", ""), os.path.basename(result[1]).replace(".fasta", ""), 1 - float(result[2]))

def get_options():
    parser = argparse.ArgumentParser(description='Create GML file from FASTA sequences using Mash distances.')
    parser.add_argument('input_fasta', help='Input FASTA file.')
    parser.add_argument('input_json', help='Input JSON file with metadata.')
    parser.add_argument('output_csv', help='Output CSV file.')
    parser.add_argument('output_gml', help='Output GML file.')
    parser.add_argument('mash_path', help='Path to mash binary.')
    parser.add_argument('--threshold', type=float, default=0.7, help='Jaccard similarity threshold (default: 0.7).')
    parser.add_argument('--threads', type=int, default=1, help='Threads for sketching (default: 1)')
    args = parser.parse_args()
    return args

def main():
    args = get_options()
    sequence_records = list(SeqIO.parse(args.input_fasta, "fasta"))
    metadata = read_metadata(args.input_json)
    # Sketch sequences using Mash
    with TemporaryDirectory() as temp_dir:
        sequence_files = Parallel(n_jobs=args.threads)(delayed(write_sequence_to_temp_file)(seq, temp_dir) for seq in tqdm(sequence_records))
        # Write sequence file paths to a text file
        sequence_files_path = os.path.join(temp_dir, "sequence_files.txt")
        with open(sequence_files_path, 'w') as f:
            for seq_file in sequence_files:
                f.write(seq_file + '\n')
        # Sketch sequences using Mash
        sys.stderr.write("Running mash sketch\n")
        sketch_file = run_mash_sketch(sequence_files_path, args.threads, args.mash_path)
        # Calculate distances
        sequence_ids = [seq.id for seq in SeqIO.parse(args.sequences, "fasta")]
        pairwise_combinations = list(combinations(sequence_ids, 2))
        distances = Parallel(n_jobs=args.jobs)(
            delayed(run_mash_dist)(seq1, seq2, sketch_file, args.mash_path)
            for seq1, seq2 in pairwise_combinations
        )
        # Filter distances
        distances = [dist for dist in distances if dist[2] <= args.threshold]

        sys.stderr.write("Writing distance GML\n")
        mash_df = pd.DataFrame(distances, columns=['seq1', 'seq2', 'distance'])
        mash_df.to_csv(args.output_csv, index=False)
        write_distance_gml(sequence_records,
                        mash_df,
                        metadata,
                        args.output_gml)
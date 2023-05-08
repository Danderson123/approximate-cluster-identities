import argparse
import json
from Bio import SeqIO

def read_mmseqs_cluster_tsv(file_path):
    cluster_assignments = {}
    clusterMapping = {}
    clusterId = 0
    with open(file_path, 'r') as f:
        for line in f:
            cluster_name, sequence_id = line.strip().split('\t')
            if not cluster_name in clusterMapping:
                clusterMapping[cluster_name] = clusterId
                clusterId += 1
            cluster_assignments[sequence_id] = clusterMapping[cluster_name]
    return cluster_assignments

def write_cluster_json(cluster_assignments, fasta_file, output_file):
    with open(output_file, 'w') as f:
        sequence_clusters = {}
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            sequence_id = seq_record.id
            if sequence_id in cluster_assignments:
                sequence_clusters[sequence_id] = cluster_assignments[sequence_id]
        json.dump(sequence_clusters, f)

def main():
    parser = argparse.ArgumentParser(description='Process MMseqs easy-cluster TSV output and FASTA file to create a JSON file with cluster assignments.')
    parser.add_argument('tsv_file', help='Input MMseqs easy-cluster TSV file.')
    parser.add_argument('fasta_file', help='Input FASTA file.')
    parser.add_argument('output_file', help='Output JSON file with cluster assignments.')
    args = parser.parse_args()

    cluster_assignments = read_mmseqs_cluster_tsv(args.tsv_file)
    write_cluster_json(cluster_assignments, args.fasta_file, args.output_file)

if __name__ == "__main__":
    main()

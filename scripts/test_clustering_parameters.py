import argparse
from Bio import SeqIO
import json
import os
import subprocess
from tqdm import tqdm

def get_options():
    parser = argparse.ArgumentParser(description='Cluster a test dataset using a clustering tool and then assess the clusterig with ACI')
    parser.add_argument('--mmseqs2', dest='mmseqs2', action='store_true', default=False,
                        help='Evaluate mmmseqs2 clusering.')
    parser.add_argument('--cdhit', dest='cdhit', action='store_true', default=False,
                        help='Evaluate cdhit clusering.')
    args = parser.parse_args()
    return args

if __name__=="__main__":

    fasta_file = "test/test_sequences.fa"
    threads = 8
    output_dir = "cluster_tool_tests"
    min_seq_id = 0.98
    temp_dir = "/tmp"
    mmseqs_bin_path = "mmseqs/bin/mmseqs"
    cdhit_bin_path = "cdhit"

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    args = get_options()

    if args.mmseqs2:
        if not os.path.exists(os.path.join(output_dir, "mmseqs2")):
            os.mkdir(os.path.join(output_dir, "mmseqs2"))
        for cov_mode in tqdm(["5", "0", "1"]):
            if not os.path.exists(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}")):
                os.mkdir(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}"))
            for c in ["0.1", "0.5","0.9"]:
                if not os.path.exists(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}")):
                    os.mkdir(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}"))
                cluster_command = " ".join([mmseqs_bin_path,
                                    "easy-cluster",
                                    fasta_file,
                                    os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}"),
                                    "/tmp",
                                    "--min-seq-id",
                                    str(min_seq_id),
                                    "--kmer-per-seq",
                                    "100",
                                    "--cov-mode",
                                    cov_mode,
                                    "-c",
                                    c,
                                    "--threads",
                                    str(threads)])
                subprocess.run(cluster_command,
                            shell=True,
                            check=True)
                metadata_command = " ".join(["python3",
                                            "scripts/get_mmseqs_cluster_assignments.py",
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster.tsv"),
                                            fasta_file,
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_assignments.json")])
                subprocess.run(metadata_command,
                            shell=True,
                            check=True)
                visualise_command = " ".join(["python3",
                                            "approximate_cluster_identities/__main__.py",
                                            #"--clusterGML",
                                            #os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_identities"),
                                            "--clusterPlot",
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_identities"),
                                            "--kmerSize",
                                            "9",
                                            "--windowSize",
                                            "50",
                                            "--threshold",
                                            str(min_seq_id),
                                            "--threads",
                                            str(threads),
                                            fasta_file,
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_assignments.json")])
                subprocess.run(visualise_command,
                            shell=True,
                            check=True)
                # visualise_command = " ".join(["python3",
                #                             "approximate_cluster_identities/__main__.py",
                #                             #"--clusterGML",
                #                             #os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_identities_relative_to_shorter"),
                #                             "--clusterPlot",
                #                             os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_identities_relative_to_shorter"),
                #                             "--kmerSize",
                #                             "9",
                #                             "--windowSize",
                #                             "50",
                #                             "--threshold",
                #                             str(min_seq_id),
                #                             "--threads",
                #                             str(threads),
                #                             "--shorter",
                #                             fasta_file,
                #                             os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_assignments.json")])
                # subprocess.run(visualise_command,
                #             shell=True,
                #             check=True)
                # rerun with cluster_mode 1
                cluster_command = " ".join([mmseqs_bin_path,
                                    "easy-cluster",
                                    fasta_file,
                                    os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1"),
                                    "/tmp",
                                    "--min-seq-id",
                                    str(min_seq_id),
                                    "--kmer-per-seq",
                                    "100",
                                    "--cov-mode",
                                    cov_mode,
                                    "-c",
                                    c,
                                    "--threads",
                                    str(threads),
                                    "--cluster-mode",
                                    "1"])
                subprocess.run(cluster_command,
                            shell=True,
                            check=True)
                metadata_command = " ".join(["python3",
                                            "scripts/get_mmseqs_cluster_assignments.py",
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster.tsv"),
                                            fasta_file,
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster_assignments.json")])
                subprocess.run(metadata_command,
                            shell=True,
                            check=True)
                visualise_command = " ".join(["python3",
                                            "approximate_cluster_identities/__main__.py",
                                            #"--clusterGML",
                                            #os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster-identities"),
                                            "--clusterPlot",
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster-identities"),
                                            "--kmerSize",
                                            "9",
                                            "--windowSize",
                                            "1",
                                            "--threshold",
                                            str(min_seq_id),
                                            "--threads",
                                            str(threads),
                                            fasta_file,
                                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster_assignments.json")])
                subprocess.run(visualise_command,
                            shell=True,
                            check=True)
                # visualise_command = " ".join(["python3",
                #                             "approximate_cluster_identities/__main__.py",
                #                             #"--clusterGML",
                #                             #os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster_identities_relative_to_shorter"),
                #                             "--clusterPlot",
                #                             os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster_identities_relative_to_shorter"),
                #                             "--kmerSize",
                #                             "9",
                #                             "--windowSize",
                #                             "1",
                #                             "--threshold",
                #                             str(min_seq_id),
                #                             "--threads",
                #                             str(threads),
                #                             "--shorter",
                #                             fasta_file,
                #                             os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_mode_1_cluster_assignments.json")])
                # subprocess.run(visualise_command,
                #             shell=True,
                #             check=True)

    if args.cdhit:
        if not os.path.exists(os.path.join(output_dir, "cdhit")):
            os.mkdir(os.path.join(output_dir, "cdhit"))
        # we have to modify the sequence headers because CDHIT cannot handle them if they're too long
        header_mapping = {}
        new_headers = []
        seqId = 0
        seq = SeqIO.parse(open(fasta_file), "fasta")
        for record in seq:
            new_headers.append(">" + str(seqId) + "\n" + str(record.seq))
            header_mapping[seqId] = record.id
            seqId += 1
        with open(os.path.join(output_dir, "cdhit", "modified_headers.fasta"), "w") as o:
            o.write("\n".join(new_headers))
        with open(os.path.join(output_dir, "cdhit", "modified_header_mappings.json"), "w") as o:
            o.write(json.dumps(header_mapping))
        for c in tqdm(["0.1", "0.5", "0.9"]):
            if not os.path.exists(os.path.join(output_dir, "cdhit", f"cov_{c}")):
                os.mkdir(os.path.join(output_dir, "cdhit", f"cov_{c}"))
            cluster_command = " ".join([cdhit_bin_path,
                                    "-i",
                                    os.path.join(output_dir, "cdhit", "modified_headers.fasta"),
                                    "-o",
                                    os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}"),
                                    "-c",
                                    str(min_seq_id),
                                    "-T",
                                    str(threads),
                                    "-s",
                                    str(c),
                                    "-n",
                                    "4",
                                    "-g",
                                    "1",
                                    "-M",
                                    "15000"])
            subprocess.run(cluster_command,
                        shell=True,
                        check=True)
            metadata_command = " ".join(["python3",
                                    "scripts/get_cdhit_cluster_assignments.py",
                                    os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}.clstr"),
                                    os.path.join(output_dir, "cdhit", "modified_headers.fasta"),
                                    os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}_cluster_assignments.json")])
            subprocess.run(metadata_command,
                            shell=True,
                            check=True)
            visualise_command = " ".join(["python3",
                                        "approximate_cluster_identities/__main__.py",
                                        #"--clusterGML",
                                        #os.path.join(output_dir, "cdhit", f"cov_{id}", f"cdhit_cov_{c}_cluster_identities"),
                                        "--clusterPlot",
                                        os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}_cluster_identities"),
                                        "--kmerSize",
                                        "9",
                                        "--windowSize",
                                        "1",
                                        "--threshold",
                                        str(min_seq_id),
                                        "--threads",
                                        str(threads),
                                        os.path.join(output_dir, "cdhit", "modified_headers.fasta"),
                                        os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}_cluster_assignments.json")])
            subprocess.run(visualise_command,
                            shell=True,
                            check=True)
            # visualise_command = " ".join(["python3",
            #                             "approximate_cluster_identities/__main__.py",
            #                             #"--clusterGML",
            #                             #os.path.join(output_dir, "cdhit", f"cov_{id}", f"cdhit_cov_{c}_cluster_identities_relative_to_shorter"),
            #                             "--clusterPlot",
            #                             os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}_cluster_identities_relative_to_shorter"),
            #                             "--kmerSize",
            #                             "9",
            #                             "--windowSize",
            #                             "1",
            #                             "--threshold",
            #                             str(min_seq_id),
            #                             "--threads",
            #                             str(threads),
            #                             "--shorter",
            #                             os.path.join(output_dir, "cdhit", "modified_headers.fasta"),
            #                             os.path.join(output_dir, "cdhit", f"cov_{c}", f"cdhit_cov_{c}_cluster_assignments.json")])
            # subprocess.run(visualise_command,
            #                 shell=True,
            #                 check=True)
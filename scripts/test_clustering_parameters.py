import os
import subprocess

fasta_file = "test/test_sequences.fa"
threads = 8
output_dir = "cluster_tool_tests"
min_seq_id = 0.9
temp_dir = "/tmp"

for o in [output_dir, os.path.join(output_dir, "mmseqs2")]:
    if not os.path.exists(o):
        os.mkdir(o)

for cov_mode in ["0", "1", "5"]:
    if not os.path.exists(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}")):
        os.mkdir(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}"))
    for c in ["0.1", "0.5", "0.9"]:
        if not os.path.exists(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}")):
            os.mkdir(os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}"))
        cluster_command = " ".join(["/home/daniel/Documents/GitHub/amira_prototype/panRG_building_tools/MMseqs2/build/bin/mmseqs",
                            "easy-cluster",
                            fasta_file,
                            os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}"),
                            "/tmp",
                            "--min-seq-id",
                            str(min_seq_id),
                            "--kmer-per-seq",
                            "1000",
                            "--cov-mode",
                            cov_mode,
                            "-c",
                            c,
                            "--threads",
                            str(threads),
                            "--cluster-reassign"])
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
                                    "--clusterGML",
                                    os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster-identities"),
                                    "--clusterPlot",
                                    os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster-identities"),
                                    "--kmerSize",
                                    "9",
                                    "--windowSize",
                                    "10",
                                    "--threshold",
                                    str(min_seq_id),
                                    "--threads",
                                    str(threads),
                                    fasta_file,
                                    os.path.join(output_dir, "mmseqs2", f"cov_mode_{cov_mode}", f"cov_{c}", f"mmseqs_cov_mode_{cov_mode}_c_{c}_cluster_assignments.json")])
        subprocess.run(visualise_command,
                    shell=True,
                    check=True)
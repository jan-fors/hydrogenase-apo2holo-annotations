"""
mmseqs easy-cluster 2026-03-06_dedup-dimers.fasta clusterRes tmp --min-seq-id 1.0 -c 0.8 --cov-mode 1
"""

import os
import subprocess

input_fasta = "2026-03-06_dedup-dimers.fasta"
output_folder = "clusterRes"
os.makedirs(output_folder, exist_ok=True)

# Run mmseqs easy-cluster
subprocess.run([
    "mmseqs", "easy-cluster",
    input_fasta,
    os.path.join(output_folder, "clusterRes_095"),
    "tmp",
    "--min-seq-id", "0.95",
    "-c", "0.8",
    "--cov-mode", "1"
])
"""
Docstring for unique_hyd_structures.03_create_unique_subset

Read the cluster results from mmseqs and create a new folder where each 
structure that is a cluster representative is copied. The output folder is called 2026-03-06_dedup-dimers_unique.
"""

import os
import shutil

input_folder = "2026-03-06-dedup_dimers"
cluster_file = "clusterRes/clusterRes_08_cluster.tsv"
output_folder = "2026-03-06_dedup-dimers_unique_08"
os.makedirs(output_folder, exist_ok=True)   
# Read cluster representatives
representatives = set()
with open(cluster_file, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            representatives.add(parts[0])  # The first column is the representative

# Copy representative structures to the output folder
for rep in representatives:
    rep_file = os.path.join(input_folder, f"{rep}")
    if os.path.exists(rep_file):
        shutil.copy(rep_file, output_folder)
    else:
        print(f"Warning: Representative file {rep_file} not found.")


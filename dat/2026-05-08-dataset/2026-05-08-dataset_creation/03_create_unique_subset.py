"""
Docstring for unique_hyd_structures.03_create_unique_subset

Read the cluster results from mmseqs and create a new folder where each 
structure that is a cluster representative is copied. The output folder is called 2026-03-06_dedup-dimers_unique.
"""

import os
import shutil

def copy_and_rename(src_path, dest_path,old_name, new_name):
	# Copy the file
	shutil.copy(src_path, dest_path)

	# Rename the copied file
	new_path = f"{dest_path}/{new_name}"
	shutil.move(f"{dest_path}/{old_name}", new_path)

input_folder = "2026-05-08-raw_structures_cleaned"
cluster_file = "result2_cluster.tsv"
output_folder = "2026-05-08-train-data_08"
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
    rep_file = os.path.join(input_folder, f"{rep}.pdb")
    if os.path.exists(rep_file):
        copy_and_rename(rep_file, output_folder,f"{rep}.pdb", f"{rep.split('_')[0]}.pdb")
    else:
        print(f"Warning: Representative file {rep_file} not found.")


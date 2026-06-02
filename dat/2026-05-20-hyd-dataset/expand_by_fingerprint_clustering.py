"""

"""

import argparse
from pathlib import Path
import pandas as pd
import csv
from sklearn.metrics import pairwise_distances
import numpy as np
from sklearn.cluster import AffinityPropagation
import shutil
import os

FACTOR=1.0

def expand(tm_cluster, feature_table, src, dest):
    """
    """
    # read foldseek results
    tm_cl = pd.read_csv(tm_cluster, sep="\t", header=None)
    tm_cl.columns = ["representatives", "samples"]

    # read features
    rows = []
    with open(feature_table, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    df = pd.DataFrame(rows) 
    rcsb = []
    for i in range(df.shape[0]):
        code = df.iloc[i]["input_path"].split("_")[0]
        rcsb.append(code)
    df["rcsb"] = rcsb

    df["input_path"] = df["input_path"].str[:-4]

    # merge
    for i in range(df.shape[0]):
        df.loc[i, "tm_repr"] = tm_cl[tm_cl["samples"] == df.loc[i, "input_path"]].iloc[0]["representatives"]

    unique_representatives = df['tm_repr'].unique() 
    print(f"unique tm_repr {unique_representatives}")

    expanded = []

    for s in unique_representatives:
        print("="*6, s, "="*6)

        # create subset
        subset = df[df["tm_repr"] == s]
        subset = subset[["input_path", "0_cl_fingerprint", "1_cl_fingerprint", "2_cl_fingerprint"]].reset_index(drop=True)

        print(f"tm_representative_cluster {s} has {subset.shape[0]} samples")

        # create input
        X = []

        for i in range(subset.shape[0]):
            cl0 = subset.loc[i, "0_cl_fingerprint"].split(",")
            cl0 = [int(x) for x in cl0]

            cl1 = subset.loc[i, "1_cl_fingerprint"].split(",")
            cl1 = [int(x) for x in cl1]

            cl2 = subset.loc[i, "2_cl_fingerprint"].split(",")
            cl2 = [int(x) for x in cl2]
            
            fs = cl0 + cl1 + cl2
            
            X.append(fs)

        sim_matrix = -pairwise_distances(X, metric="euclidean")

        print(f"median: {np.median(sim_matrix):.2f}")
        print(f"min:    {np.min(sim_matrix):.2f}")

        # start at min to get coarser clusters
        af = AffinityPropagation(
            preference=FACTOR*np.min(sim_matrix),
            damping=0.9,                    
            max_iter=500,
            convergence_iter=50
        ).fit(X)

        print(f"n_clusters: {len(af.cluster_centers_indices_)}")
        
        for i in af.cluster_centers_indices_:
            expanded.append(subset.iloc[i]["input_path"])

    print(f"Sample Size after Clustering: {len(expanded)}")

    if not os.path.exists(dest):
        os.makedirs(dest, exist_ok=True)

    for e in expanded:
        shutil.copy2(src / Path(e + ".pdb"), dest / Path(e.split("_")[0] + ".pdb"))
        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tm_cluster", type=Path, help="Path to the foldseek result file tsv.")
    parser.add_argument("feature_table", type=Path, help="Path to the extracted feature table.")
    parser.add_argument("src_dir", type=Path, help="Directory with the structures")
    parser.add_argument("dest_dir", type=Path, help="Destination directory")
    
    args = parser.parse_args()

    expand(args.tm_cluster, args.feature_table, args.src_dir, args.dest_dir)
"""
Script that creates the subsets for a k-fold cv
"""
from pathlib import Path
import os
import argparse
import random
import shutil
import subprocess
import pandas as pd
import time

def _read_test_data(test_data : Path):
    """
    """
    if not os.path.exists(test_data):
        raise ValueError(f"{test_data} does not exist")
    
    files = os.listdir(test_data)
    
    return [test_data/Path(x) for x in files]

def _create_randomized_subsets(test_data : Path, k : int):
    """
    """
    files = _read_test_data(test_data)
    
    random.shuffle(files)

    subsets = {}

    if k == -1:
        k = len(files)
    elif len(files) < k:
        k = len(files)

    for i in range(len(files)):
        if i%k not in subsets.keys():
            subsets[i%k] = []
        subsets[i%k].append(files[i])

    avg_subset_size = len(files) / k

    print("Avg. subset size:", avg_subset_size)

    return subsets

def _search_test_files_against_train_data(test_dir : Path, train_data : Path, min_seq_id : float, output_file : Path, tmp : Path):
    """
    """
    # run against foldseek easy-multimersearch benchmark_dir/0/test dat/2026-05-26-cof-dataset/seq_repr/ output delete --min-seq-id 0.95
    cmd = ["foldseek", "easy-multimersearch", str(test_dir), str(train_data), str(output_file), str(tmp), "--min-seq-id", str(min_seq_id)]

    # run
    subprocess.run(cmd, stdout=subprocess.DEVNULL)

    time.sleep(5)

    # parser result
    try:
        df = pd.read_csv(str(output_file) + "_report", sep="\t", header=None)
        hits = df.iloc[:, 1].to_list()

        for i in range(len(hits)):
            if "_" in hits[i]:
                hits[i] = hits[i].split("_")[0]

        hits = [x + '.pdb' for x in hits]

        return hits
    except Exception as e:
        return []

def create_subsets(test_data_dir : Path, training_data_dir : Path, k : int, min_seq_id : float, out : Path, tmp : Path):
    """
    """
    # read test_data and create k subsets
    subsets = _create_randomized_subsets(test_data_dir, k)

    for subset in subsets.keys():
        print("="*20,"SUBSET", subset, "="*20)
        # create subset dir
        subset_dir = out / Path(str(subset))
        os.makedirs(subset_dir, exist_ok=True)

        # create test dir
        test_dir = subset_dir / Path('test')
        os.makedirs(test_dir, exist_ok=True)
        counter = 0
        for file in subsets[subset]:
            shutil.copy2(file, test_dir/file.name)
            counter += 1

        print(f"Test-Set contains {counter} files.")
        foldseek_dir = subset_dir / Path('foldseek')
        os.makedirs(foldseek_dir, exist_ok=True)
        result_file = foldseek_dir / Path('res')

        # search test files against train_data 
        blacklisted = _search_test_files_against_train_data(test_dir, training_data_dir, min_seq_id, result_file, tmp)
        print(f"Training-Set contains {len(blacklisted)} files that are above min-seq-id {min_seq_id}")

        # create train dir
        train_dir = subset_dir / Path('train')
        os.makedirs(train_dir, exist_ok=True)

        counter = 0
        for file in os.listdir(training_data_dir):
            if file in blacklisted:
                pass
            else:
                counter += 1
                shutil.copy2(training_data_dir/Path(file), train_dir/Path(file))

        print(f"Training-Set contains {counter} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("test_data", type=Path, help="Test data src dir")
    parser.add_argument("training_data", type=Path, help="Training data src dir")
    parser.add_argument("--k", type=int, default=5, help="k-fold cv")
    parser.add_argument("--min_seq_id", type=float, default=1.0, help="min seq id for deduplicating dataset")
    parser.add_argument("--tmp", type=Path, help="Path to the tmp folder")
    parser.add_argument("--out", type=Path, help="Output directory", default="out")

    args = parser.parse_args()


    if not os.path.exists(args.test_data):
        raise ValueError("Provide valid test data path")

    if not os.path.exists(args.training_data):
        raise ValueError("Provide valid training data path")

    if not os.path.exists(args.out):
        os.makedirs(args.out, exist_ok=True)

    if not os.path.exists(args.tmp):
        os.makedirs(args.tmp)

    create_subsets(
        test_data_dir=args.test_data,
        training_data_dir=args.training_data,
        k=args.k,
        min_seq_id=args.min_seq_id,
        out=args.out,
        tmp=args.tmp
        )
"""
Benchmarking

1. Split the test structures into k subsets -> LOOCV (k=-1)?
2. search test structures against training dataset and remove copies or direct homologs
3. train and test
"""
import os
from pathlib import Path
import argparse
import random
import shutil
import subprocess
import pandas as pd
from src.build import build
from src.main import main

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

    return subsets

def _search_test_files_against_train_data(test_dir : Path, train_data : Path, min_seq_id : float, output_file : Path, tmp : Path):
    """
    """
    # run against foldseek easy-multimersearch benchmark_dir/0/test dat/2026-05-26-cof-dataset/seq_repr/ output delete --min-seq-id 0.95
    cmd = ["foldseek", "easy-multimersearch", str(test_dir), str(train_data), str(output_file), str(tmp), "--min-seq-id", str(min_seq_id)]

    # run
    subprocess.run(cmd)

    # parser result
    df = pd.read_csv(str(output_file) + "_report", sep="\t", header=None)
    hits = df.iloc[:, 1].to_list()

    for i in range(len(hits)):
        if "_" in hits[i]:
            hits[i] = hits[i].split("_")[0]

    hits = [x + '.pdb' for x in hits]

    return hits

def benchmark(test_data : Path, train_data : Path, k : int, benchmark_dir : Path, min_seq_id : float, prediction_only : bool, result_dir : Path = Path('out')):
    """

    TODO prediction_only function


    1. split the test data files into k subsets.
    2. for each subset:
        a. create a directory in benchmark dir
        b. copy the test_set files to benchmark_dir/k/test
        c. search the test_set files against the train_data
        d. parse results
        e. copy all that are not hits into benchmark_dir/k/train
        f. build database with those files
        g. predict test set files
        (h. compare raw hydrogenase-feature-extracts from test and predset) <- outside of benchmark file
    """
    if prediction_only: # assumes all folders already exist
        tmp_folder = benchmark_dir / Path('tmp')

        for subset_folder in os.listdir(benchmark_dir):
            if subset_folder == "tmp":
                continue
            subset_dir = benchmark_dir / Path(str(subset_folder))
            test_dir = subset_dir / Path('test')

            # load database
            db_dir = subset_dir / Path('db')
            chain_dict = db_dir / Path('single_chains')
            structure_db = db_dir / Path('structureDB')
            structure_db = structure_db / Path('structureDB')
            fingerprint_db = db_dir / Path('fingerprintDB')

            pred_dir = subset_dir / result_dir
           
            for test_structure in os.listdir(test_dir):
                try:
                    print("Predicting test structure", test_structure)
                    test_structure_path = test_dir / Path(test_structure)
                    output = pred_dir / Path(test_structure.split(".")[0])
                    os.makedirs(output, exist_ok=True)
                    main(input_path=test_structure_path,
                        out=output,
                        tmp=tmp_folder,
                        boltz=False,
                        plot=True, 
                        structure_db_path=structure_db,
                        fingerprint_db_path=fingerprint_db)
                except FileNotFoundError as e:
                    print(e)
    else:
        # create benchmark dir
        os.makedirs(benchmark_dir, exist_ok=True)

        # create tmp folder inside of benchmark
        tmp_folder = benchmark_dir / Path('tmp')
        os.makedirs(tmp_folder, exist_ok=True)

        # read test_data and create k subsets
        subsets = _create_randomized_subsets(test_data, k)

        for subset in subsets.keys():
            
            # create subset dir
            subset_dir = benchmark_dir / Path(str(subset))
            os.makedirs(subset_dir, exist_ok=True)

            # create test dir
            test_dir = subset_dir / Path('test')
            os.makedirs(test_dir, exist_ok=True)
            for file in subsets[subset]:
                shutil.copy2(file, test_dir/file.name)

            foldseek_dir = subset_dir / Path('foldseek')
            os.makedirs(foldseek_dir, exist_ok=True)
            result_file = foldseek_dir / Path('res')

            # search test files against train_data 
            blacklisted = _search_test_files_against_train_data(test_dir, train_data, min_seq_id, result_file, tmp_folder)

            # create train dir
            train_dir = subset_dir / Path('train')
            os.makedirs(train_dir, exist_ok=True)

            for file in os.listdir(train_data):
                if file in blacklisted:
                    #print("blacklisted")
                    pass
                else:
                    shutil.copy2(train_data/Path(file), train_dir/Path(file))

            # create database
            db_dir = subset_dir / Path('db')
            os.makedirs(db_dir, exist_ok=True)
            chain_dict = db_dir / Path('single_chains')
            structure_db = db_dir / Path('structureDB')
            os.makedirs(structure_db, exist_ok=True)
            structure_db = structure_db / Path('structureDB')
            fingerprint_db = db_dir / Path('fingerprintDB')
            os.makedirs(fingerprint_db, exist_ok=True)

            build(train_dir, chain_dict, structure_db, fingerprint_db)

            # predict
            pred_dir = subset_dir / Path('out')
            os.makedirs(pred_dir, exist_ok=True)

            for test_structure in os.listdir(test_dir):
                try:
                    print("Predicting test structure", test_structure)
                    test_structure_path = test_dir / Path(test_structure)
                    output = pred_dir / Path(test_structure.split(".")[0])
                    os.makedirs(output, exist_ok=True)
                    main(input_path=test_structure_path,
                        out=output,
                        tmp=tmp_folder,
                        boltz=False,
                        plot=True, 
                        structure_db_path=structure_db,
                        fingerprint_db_path=fingerprint_db)
                except (FileNotFoundError, subprocess.CalledProcessError) as e:
                    print(e)
                
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("test_data", type=Path, help="Test data folder")
    parser.add_argument("train_data", type=Path, help="Training data dir")
    parser.add_argument("--k", type=int, default=5, help="Split Test data into k subsets.")
    parser.add_argument("--benchmark_dir", type=Path, default="benchmark_dir", help="directory where the benchmarking subsets are created etc.")
    parser.add_argument("--min_seq_id", type=float, default=0.95, help="Minimum Sequence Identity for the foldseek run")
    parser.add_argument("--prediction_only", action="store_true")
    parser.add_argument("--result_dir_name", type=Path, default=Path('out'))

    args = parser.parse_args()

    benchmark(args.test_data, args.train_data, args.k, args.benchmark_dir, args.min_seq_id, args.prediction_only, args.result_dir_name)
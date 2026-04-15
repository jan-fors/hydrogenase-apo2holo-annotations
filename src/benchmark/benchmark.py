"""
predict each of the test cases and compare to gt
then calculate a score
"""
from pathlib import Path
import os
from src.main import main
from src.build import build
import random
import shutil
from datetime import datetime
from src.benchmark.clean_test_structures import clean_test_structures

########################### PARAMS###########################
RAW_STRUCTURE_DIR = Path("dat/2026-03-06-dedup_dimers")
BENCHMARKING_DIR = Path("benchmarking")
TARGET_DIR = Path("dat/benchmarking_subsets")
os.makedirs(BENCHMARKING_DIR, exist_ok=True)

K = 10

CREATE_SUBSETS = False
CREATE_TEST_FOLDERS = False
CREATE_DATABASES = False
PREDICT = True
MEASURE = False

# 0 create groundtruth dataset

# 1. create subsets
if CREATE_SUBSETS:
    print(f"[{datetime.now()}] Creating Benchmarking Subsets...")

    #receive all files
    files = [f for f in RAW_STRUCTURE_DIR.iterdir() if f.is_file()]
    print(f"[{datetime.now()}] Read {len(files)} files")

    # zufällig mischen
    random.shuffle(files)

    # Unterordner erstellen
    set_dirs = []
    for i in range(K):
        d = TARGET_DIR / f"set_{i+1}"
        d.mkdir(parents=True, exist_ok=True)
        set_dirs.append(d)

    # Dateien möglichst gleichmäßig verteilen
    for i, file in enumerate(files):
        dest = set_dirs[i % K] / file.name
        shutil.copy2(file, dest)

    for set_dir in os.listdir(TARGET_DIR):
        set_dir_path = Path(os.path.join(TARGET_DIR, set_dir))
        files = [f for f in set_dir_path.iterdir() if f.is_file()]
        print(f"[{datetime.now()}] {set_dir} {len(files)}")

    print(f"[{datetime.now()}] Creating Benchmarking Subsets done")

# 2. create test folders
if CREATE_TEST_FOLDERS:
    print(f"[{datetime.now()}] Creating Test Folders ...")

    for set in os.listdir(TARGET_DIR):
        testing_set_directory_path = os.path.join(BENCHMARKING_DIR, set) # benchmarking/setX/...
        os.makedirs(testing_set_directory_path, exist_ok=True)

        db_folder_path = os.path.join(testing_set_directory_path, "db") # benchmarking/setX/db
        os.makedirs(db_folder_path, exist_ok=True)
        test_set_path = Path(os.path.join(testing_set_directory_path, "test")) # benchmarking/setX/test
        os.makedirs(test_set_path, exist_ok=True)
        raw_data_path = Path(os.path.join(testing_set_directory_path, "raw_data"))
        os.makedirs(raw_data_path, exist_ok=True)

        # copy files to raw_data_path
        for in_set in os.listdir(TARGET_DIR):
            if set != in_set:
                files = [f for f in Path(os.path.join(TARGET_DIR, in_set)).iterdir() if f.is_file()]
                for i, file in enumerate(files):
                    dest = raw_data_path / file.name
                    shutil.copy2(file, dest)

        # copy files to test
        files = [f for f in Path(os.path.join(TARGET_DIR, set)).iterdir() if f.is_file()]
        for i, file in enumerate(files):
            dest = test_set_path / file.name
            shutil.copy2(file, dest)

        # clean test files
        clean_test_structures(test_set_path)

    print(f"[{datetime.now()}] Creating Test Folders done")

# 3. build databases
if CREATE_DATABASES:
    print(f"[{datetime.now()}] Building Databases ...")

    # iterate over each set
    for set in os.listdir(BENCHMARKING_DIR):
        print(f"Creating database for subset {set}")
        input_structure_dir = os.path.join(BENCHMARKING_DIR, set,  "raw_data")
        output_structure_dir = os.path.join(BENCHMARKING_DIR,set, "single_chains")
        db_dir = os.path.join(BENCHMARKING_DIR,set, "db")
        structure_db = os.path.join(db_dir, "structureDB")
        os.makedirs(structure_db, exist_ok=True)
        structure_db = os.path.join(structure_db, "structureDB")
        
        fingerprint_db = os.path.join(db_dir, "fingerprintDB_6")
        os.makedirs(fingerprint_db, exist_ok=True)

        build(input_structure_dir, output_structure_dir, structure_db, fingerprint_db)

    print(f"[{datetime.now()}] Building Databases done")

# 4. run
if PREDICT:
    print(f"[{datetime.now()}] Annotate Test Sets ...")

    for set in os.listdir(BENCHMARKING_DIR):
        structure_db_path = os.path.join(BENCHMARKING_DIR, set, "db", "structureDB", "structureDB")
        fingerprint_db = os.path.join(BENCHMARKING_DIR, set, "db", "fingerprintDB_6")

        output_dir = os.path.join(BENCHMARKING_DIR, set, "logreg_6")
        os.makedirs(output_dir, exist_ok=True)

        test_set = os.path.join(BENCHMARKING_DIR, set, "test")
        for structure in os.listdir(test_set):
            structure_path = os.path.join(test_set, structure)
            structure_name = structure.split(".")[0]

            try:
                main(input_path=structure_path, output=structure_name, output_dir=output_dir,structure_db_path=structure_db_path, fingerprint_db_path=fingerprint_db, tmp="tmp", boltz=False)
            except:
                continue

        

    print(f"[{datetime.now()}] Annotation done")

# 5. measure
if MEASURE:
    print(f"[{datetime.now()}] Measure Accuracy ...")



    print(f"[{datetime.now()}] Measurement done")

# 6. result
print("############# RESULTS #############")

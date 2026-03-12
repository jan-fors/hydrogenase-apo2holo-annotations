"""
using a k-fold cv 
create 10 structure databases and 10 fingerprint databases (or fingerprint model depending on approach here)
"""
from pathlib import Path
import os
import shutil
from src.database.build_structure_db import build_structure_db
from src.database.build_fingerprint_db import build_fingerprint_db

BENCHMARKING_DB_PARENT = Path("db/benchmarking_dbs")
DATASET_DIR = Path("dat/benchmarking_subsets")

# build_structure_db
for set in os.listdir(DATASET_DIR):
    out_path = os.path.join(BENCHMARKING_DB_PARENT, set)
    os.makedirs(out_path, exist_ok=True)

    data_dir = os.path.join(out_path, "database_files")
    os.makedirs(data_dir, exist_ok=True)
    
    data2_dir = os.path.join(out_path, "chain_files")
    os.makedirs(data_dir, exist_ok=True)

    test_data = os.path.join(out_path, "test")
    os.makedirs(test_data, exist_ok=True)

    structure_db_path = os.path.join(out_path, "structureDB")

    for set2 in os.listdir(DATASET_DIR):
        if set != set2:
            for file in os.listdir(os.path.join(DATASET_DIR, set2)):
                shutil.copy(os.path.join(DATASET_DIR, set2, file), os.path.join(data_dir, file))
        else:

            for file in os.listdir(os.path.join(DATASET_DIR, set2)):
                shutil.copy(os.path.join(DATASET_DIR, set2, file), os.path.join(test_data, file))

    build_structure_db(data_dir, data2_dir, structure_db_path)
    build_fingerprint_db(data_dir, out_path)

    
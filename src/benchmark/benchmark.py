"""
predict each of the test cases and compare to gt
then calculate a score
"""
from pathlib import Path
import os
from src.main import main

BENCHMARK_DIR = Path("/home/solar/Documents/hiwi/deephyds/hydrogenase-apo2holo-annotations/db/benchmarking_dbs")

for set in os.listdir(BENCHMARK_DIR):
    test_set = os.path.join(BENCHMARK_DIR, set, "test")
    structure_db = os.path.join(BENCHMARK_DIR, set, "structureDB")
    fingerprint_db = os.path.join(BENCHMARK_DIR, set, "fingerprint_db.tsv")
    print(set)
    for apo in os.listdir(test_set):
        structure_path = os.path.join(test_set, apo)
        outputs = os.path.join(BENCHMARK_DIR, set, "out")

        print(apo)
        main(structure_path, apo, outputs, os.path.join(outputs, "tmp"), False, structure_db, fingerprint_db)
"""
python -m src.benchmark.complete_benchmark.deduplicate_fingeprint_db benchmarking/complete_benchmark/subsets db_R5.0.tsv
"""

import argparse
import os
import pandas as pd
from pathlib import Path
from apo2holo.database.utils.deduplicator import deduplicate


def deduplicate_fingerprint_db(
    directory_path: Path, db_name: str, output_name_prefix: str = "dd_"
) -> str:
    """
    Iterates over all subsets and deduplicates the database by looking at the fingerprints
    Returns the name of the new database
    """
    for subset in os.listdir(directory_path):
        # define db path
        db_path = (
            directory_path
            / Path(subset)
            / Path("db")
            / Path("fingerprintDB")
            / Path(db_name)
        )

        # read db
        if str(db_path).endswith(".gz"):
            df = pd.read_csv(db_path, sep="\t", compression="gzip")
        else:
            df = pd.read_csv(db_path, sep="\t")

        # deduplicate
        new_df = deduplicate(df)

        # define output path
        output_file_name = output_name_prefix + db_name
        out_path = (
            directory_path
            / Path(subset)
            / Path("db")
            / Path("fingerprintDB")
            / Path(output_file_name)
        )

        # save
        if str(out_path).endswith(".gz"):
            new_df.to_csv(out_path, index=False, sep="\t", compression="gzip")
        else:
            new_df.to_csv(out_path, index=False, sep="\t")

    return output_file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("db_name", type=str)

    args = parser.parse_args()

    deduplicate_fingerprint_db(args.dir, args.db_name)

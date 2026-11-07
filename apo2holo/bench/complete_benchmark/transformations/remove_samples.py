"""
Removes samples matching different criteria
E.g. Remove 4Fe4S where aug_dist > 0.0
"""
from pathlib import Path
import os
import pandas as pd
import operator


OPS = {"==": operator.eq, "!=": operator.ne, ">": operator.gt,
       ">=": operator.ge, "<": operator.lt, "<=": operator.le}


def remove_samples_from_db(directory_path : Path, config: dict, db_name : str, output_name_prefix : str = "rm_") -> str:
    """
    Iterates over all fingerprint tsvs in all subsets and applies defined filtering rules
    Returns name of transformed databaase.
    """
    for subset in os.listdir(directory_path):
        # define db path
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
    
        # read db
        if str(db_path).endswith(".gz"):
            df = pd.read_csv(db_path, sep="\t", compression="gzip")
        else:
            df = pd.read_csv(db_path, sep="\t")

        # apply filters
        for where in config["where"]:
            new_df = drop_matching(df, where)

        # define output
        output_name = output_name_prefix+db_name
        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(output_name)

        # save
        if str(out_path).endswith(".gz"):
            new_df.to_csv(out_path, index=False, sep="\t", compression="gzip")
        else:
            new_df.to_csv(out_path, index = False, sep = "\t")

    return output_name

def drop_matching(df, conditions):
    mask = pd.Series(True, index=df.index)
    for col, (op, val) in conditions.items():
        mask &= OPS[op](df[col], val)
    return df[~mask]


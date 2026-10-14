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
    """
    for subset in os.listdir(directory_path):
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
    
        # read db
        df = pd.read_csv(db_path, sep="\t")


        for where in config["where"]:

            new_df = drop_matching(df, where)


        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(output_name_prefix+db_name)
        new_df.to_csv(out_path, index = False, sep = "\t")

    return output_name_prefix+db_name

def drop_matching(df, conditions):
    mask = pd.Series(True, index=df.index)
    for col, (op, val) in conditions.items():
        mask &= OPS[op](df[col], val)
    return df[~mask]


"""
"""
from pathlib import Path
import pandas as pd
from src.utils.constants import (
    FOLDSEEK_OUT_FORMAT
)

def parse_msearch_output(res_file : Path) -> pd.DataFrame:
    """
    """
    columns = FOLDSEEK_OUT_FORMAT.split(",")
    return pd.read_csv(res_file, sep="\t", header=None, names=columns)
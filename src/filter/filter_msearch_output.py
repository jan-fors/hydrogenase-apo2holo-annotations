import pandas as pd
from src.utils.constants import (
    FIDENT_THRESHOLD,
    BITS_THRESHOLD
)

def filter_msearch_output(df : pd.DataFrame, fident_threshold : float, bits_threshold : float) -> pd.DataFrame:
    """
    """
    subset = df.copy()
    subset = subset[subset["fident"] > fident_threshold]
    subset = subset[subset["bits"] > bits_threshold]
    
    return subset
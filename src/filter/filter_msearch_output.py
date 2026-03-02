import pandas as pd
from src.utils.constants import (
    FIDENT_THRESHOLD,
    BITS_THRESHOLD
)

def filter_msearch_output(df : pd.DataFrame) -> pd.DataFrame:
    """
    """
    subset = df.copy()
    subset = subset[subset["fident"] > FIDENT_THRESHOLD]
    subset = subset[subset["bits"] > BITS_THRESHOLD]
    
    return subset
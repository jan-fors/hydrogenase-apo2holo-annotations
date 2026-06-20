import pandas as pd

def deduplicate(df : pd.DataFrame) -> pd.DataFrame:
    cols = ['ALA', 'ARG', 'ASN', 'ASP',
       'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
       'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    df["fingerprint"] = df[cols].apply(lambda row:"_".join(row.values.astype(str)), axis=1)
    subset = df.drop_duplicates(subset=["fingerprint"], keep='first')
    return subset
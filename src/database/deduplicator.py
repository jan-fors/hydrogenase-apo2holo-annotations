import pandas as pd

def deduplicate(df : pd.DataFrame) -> pd.DataFrame:
   cols = list(df.columns)

   if "res_name" in cols:
      cols.remove("res_name")

   if "id" in cols:
      cols.remove("id")

   if "formula" in cols:
      cols.remove("formula")

   if "smiles" in cols:
      cols.remove("smiles")



    
   df["fingerprint"] = df[cols].apply(lambda row:"_".join(row.values.astype(str)), axis=1)
   subset = df.drop_duplicates(subset=["fingerprint"], keep='first')
   return subset
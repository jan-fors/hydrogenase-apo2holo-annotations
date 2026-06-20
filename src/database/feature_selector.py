import pandas as pd

def feature_selector(data : pd.DataFrame, columns : list):
    return data[columns]
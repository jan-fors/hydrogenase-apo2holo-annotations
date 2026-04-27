import pandas as pd
import os
from src.utils.constants import (
    AA_ORDER
)
from collections import Counter
from pathlib import Path
import numpy as np

from src.io.printl import printl

class AbsolutSearchEngine:
    def __init__(self):
        pass

    def load_file(self, table_path : Path):
        """
        """
        if os.path.exists(table_path):
            printl(f"Loading table from {table_path}")
            self.table = pd.read_csv(table_path, sep="\t")
            self.table = self.table.set_index(AA_ORDER)
            return self
        else:
            raise ValueError(f"{model_path} does not exist.")

    def load_dataframe(self, df : pd.DataFrame):
        """
        """
        self.table = df

    def predict(self, F : np.array) -> list:
        """
        F can have multiple items
        """
        results = []
        for row in F:
            key = tuple(row)
            if key in self.table.index:
                results.extend(self.table.loc[key]["formula"])

        # parse results into probas
        counts = Counter(results)

        total = len(results)

        proba = {k: v for k, v in counts.items()}

        return proba

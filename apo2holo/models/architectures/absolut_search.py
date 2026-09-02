import pandas as pd
import os
from apo2holo.utils.constants import (
    STANDARD_AMINO_ACIDS
)
from collections import Counter
from pathlib import Path
import numpy as np

from apo2holo.io.writers.printl import printl

class AbsolutSearchEngine:
    def __init__(self):
        pass

    def load_file(self, table_path : Path):
        """
        """
        if os.path.exists(table_path):
            printl(f"Loading table from {table_path}")
            self.table = pd.read_csv(table_path, sep="\t")
            self.table = self.table.set_index(STANDARD_AMINO_ACIDS)
            return self
        else:
            raise ValueError(f"{table_path} does not exist.")

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

from sklearn.neural_network import MLPClassifier
import os
import numpy as np
from src.utils.constants import (
    MLP_ACTIVATION,
    MLP_SOLVER,
    MLP_ALPHA,
    MLP_MAX_ITER,
    MLP_BATCH_SIZE,
    MLP_BETA_1,
    MLP_BETA_2,
    MLP_EARLY_STOPPING,
    MLP_HIDDEN_LAYERS,
    MLP_LEARNING_RATE,
    MLP_LEARNING_RATE_INIT,
    MLP_MOMENTUM,
    MLP_TOL,
    MLP_VALIDATION_FRACTION
)

from pathlib import Path
import pickle
import argparse
import pandas as pd
from src.io.printl import printl

class CustomMLPClassifier:
    def __init__(self):
        pass

    def load(self, model_path : Path):
        """
        """
        if os.path.exists(model_path):
            printl(f"Loading model from {model_path}")
            with open(model_path, "rb") as f:
                self.model = pickle.load(f)
            
            return self
        else:
            raise ValueError(f"{model_path} does not exist.")

    def predict(self, F : np.array) -> list:
        """
        F can have multiple items
        """
        
        proba = self.model.predict_proba(F)

        class_names = self.model.classes_

        output = [
            {cls: p for cls, p in zip(class_names, row)}
            for row in proba
        ]

        return output

    def train(self, X, Y):
        """
        """
        self.model = MLPClassifier(
                activation=MLP_ACTIVATION,
                alpha=MLP_ALPHA,
                batch_size=MLP_BATCH_SIZE,
                beta_1=MLP_BETA_1,
                beta_2=MLP_BETA_2,
                early_stopping=MLP_EARLY_STOPPING,
                hidden_layer_sizes=MLP_HIDDEN_LAYERS,
                learning_rate=MLP_LEARNING_RATE,
                learning_rate_init=MLP_LEARNING_RATE_INIT,
                max_iter=MLP_MAX_ITER,
                momentum=MLP_MOMENTUM,
                solver=MLP_SOLVER,
                tol=MLP_TOL,
                validation_fraction=MLP_VALIDATION_FRACTION
            )
        self.model.fit(X, Y)
    
    def get_model(self):
        """
        """
        return self.model

    def save(self, path : Path):
        """
        
        """
        with open(path, "wb") as f:
            pickle.dump(self.model, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", type = Path, help="path to the fingerprint tsv file")
    parser.add_argument("name", type=str, help="outputname")
    parser.add_argument("-o", type=Path, default=Path("."))
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    X = df.drop(columns=['formula', 'id', 'smiles', 'res_name'])
    X_np = X.to_numpy()

    Y = df["formula"]
    model = CustomMLPClassifier()
    model.train(X,Y)

    model.save(args.o / Path(args.name))
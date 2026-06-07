
from pathlib import Path
import os
from src.io.printl import printl
from sklearn.linear_model import LogisticRegression
import pickle
import numpy as np
from src.utils.constants import (
    LOGREG_C,
    LOGREG_CLASS_WEIGHT,
    LOGREG_L1_RATIO,
    LOGREG_MAX_ITER,
    LOGREG_SOLVER,
    LOGREG_FIT_INTERCEPT,
    LOGREG_TOL
)

class LogisticRegressionModel:
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
        self.model = LogisticRegression(
            C=LOGREG_C,
            class_weight=LOGREG_CLASS_WEIGHT,
            fit_intercept=LOGREG_FIT_INTERCEPT,
            l1_ratio=LOGREG_L1_RATIO,
            max_iter=LOGREG_MAX_ITER,
            solver=LOGREG_SOLVER,
            tol=LOGREG_TOL
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
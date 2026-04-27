from sklearn.svm import SVC
from sklearn import svm
import os
import numpy as np
from src.utils.constants import (
    SVC_MAX_ITER,
    SVC_KERNEL,
    SVC_DEGREE,
    SVC_GAMMA,
    SVC_SHRINKING,
    SVC_TOLERANCE,
    SVC_C
)
from pathlib import Path
import pickle

from src.io.printl import printl

class SupportVectorMachine:
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
        self.model = svm.SVC(
            probability=True,
            C=SVC_C,
            degree=SVC_DEGREE,
            max_iter=SVC_MAX_ITER,
            gamma=SVC_GAMMA,
            kernel=SVC_KERNEL,
            tol=SVC_TOLERANCE,
            shrinking=SVC_SHRINKING
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
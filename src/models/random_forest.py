import os
from pathlib import Path
import pickle
from src.io.printl import printl
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from src.utils.constants import (
    RANF_BOOTSTRAP,
    RANF_CLASS_WEIGHT,
    RANF_MAX_DEPTH,
    RANF_MAX_FEATURES,
    RANF_MAX_SAMPLES,
    RANF_MIN_SAMPLE_SPLIT,
    RANF_MIN_SAMPLES_LEAF,
    RANF_N_ESTIMATORS
)

class CustomRandomForestClassifier:
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
    
    def get_model(self):
        """
        """
        return self.model

    def save(self, path : Path):
        """
        
        """
        with open(path, "wb") as f:
            pickle.dump(self.model, f)

    
    def train(self, X, Y):
        """
        """
        self.model = RandomForestClassifier(
            bootstrap=RANF_BOOTSTRAP,
            class_weight=RANF_CLASS_WEIGHT,
            max_depth=RANF_MAX_DEPTH,
            max_features=RANF_MAX_FEATURES,
            max_samples=RANF_MAX_SAMPLES,
            min_samples_leaf=RANF_MIN_SAMPLES_LEAF,
            min_samples_split=RANF_MIN_SAMPLE_SPLIT,
            n_estimators=RANF_N_ESTIMATORS
            )
        self.model.fit(X, Y)
    
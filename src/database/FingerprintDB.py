import os
from pathlib import Path
from src.utils.geometric.point_distance import point_distance
from typing import Literal
import pickle
import numpy as np
import pandas as pd
from src.io.writers.printl import printl
from src.utils.constants import (
    FINGERPRINT_DB,
    MODEL,
    STANDARD_AMINO_ACIDS,
    FINGERPRINT_RADIUS,
    MIN_DIST_ARTIFICAL_SAMPLES,
    AMOUNT_ARTIFICAL_SAMPLES,
    N_AUGMENTATIONS,
    AUGMENTATION_RADIUS
)
from sklearn.linear_model import LogisticRegression
from collections import Counter
from src.filter.apply_blacklist import apply_blacklist_build
from src.filter.apply_whitelist import apply_whitelist_build
from src.utils.geometric.calculate_geometric_centers import calculate_geometric_centers

from src.fingerprint.create_fingerprint import create_fingerprint
import uuid
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings

warnings.simplefilter("ignore", PDBConstructionWarning)
warnings.filterwarnings("ignore", message="X does not have valid feature names")


class FingerprintDB:
    def __init__(self, db_directory: Path = None):
        pass

    def as_search(
        self,
        F: dict,
        search_type: Literal[
            "logreg", "absolut", "svm", "mlp", "nn", "randforest", None, "sum", "mc"
        ],
        model: Path = None,
    ) -> list:
        """
        Return:
            TODO return the same type
        """
        if search_type == "sum" or search_type == "mc" or search_type == None:
            ################ TODO REMOVE AFTER BENCHMARKING ################
            # load model
            if os.path.exists(model):
                printl(f"Loading model from {model}")
                with open(model, "rb") as f:
                    smodel = pickle.load(f)
            else:
                raise ValueError(f"{model} does not exist.")
            # predict
            proba = smodel.predict_proba(F)

            class_names = smodel.classes_

            output = [{cls: p for cls, p in zip(class_names, row)} for row in proba]
            return output
            ################ TODO REMOVE AFTER BENCHMARKING ################
        else:
            if search_type == "logreg":
                return self.as_logistic_regression_model.predict(F)
            elif search_type == "absolut":
                return self.as_absolut_search_engine.predict(F)
            elif search_type == "svm":
                return self.as_svm_model.predict(F)
            elif search_type == "mlp":
                return self.as_mlp_classifier_model.predict(F)
            elif search_type == "randforest":
                return self.as_random_forest_model.predict(F)
            # elif search_type == "nn":
            #     return self._nn_predict(F)
            else:
                printl("Use valid earch option. [logreg, absolut]")
                return []

    def fes_type_search(
        self,
        F: dict,
        search_type: Literal[
            "logreg", "absolut", "svm", "mlp", "nn", "randforest", None, "sum", "mc"
        ],
        model: Path = None,
    ) -> list:
        """
        Return:
            TODO return the same type
        """
        if search_type == "sum" or search_type == "mc" or search_type == None:
            ################ TODO REMOVE AFTER BENCHMARKING ################
            # load model
            if os.path.exists(model):
                printl(f"Loading model from {model}")
                with open(model, "rb") as f:
                    smodel = pickle.load(f)
            else:
                raise ValueError(f"{model} does not exist.")
            # predict
            proba = smodel.predict_proba(F)

            class_names = smodel.classes_

            output = [{cls: p for cls, p in zip(class_names, row)} for row in proba]
            return output
            ################ TODO REMOVE AFTER BENCHMARKING ################
        else:
            if search_type == "logreg":
                return self.fes_logistic_regression_model.predict(F)
            elif search_type == "absolut":
                return self.fes_absolut_search_engine.predict(F)
            elif search_type == "svm":
                return self.fes_svm_model.predict(F)
            elif search_type == "mlp":
                return self.fes_type_mlp_classifier_model.predict(F)
            elif search_type == "randforest":
                return self.fes_random_forest_model.predict(F)
            # elif search_type == "nn":
            #     return self._nn_predict(F)
            else:
                printl("Use valid earch option. [logreg, absolut]")
                return []

    def fes_pocket_search(
        self,
        F: dict,
        search_type: Literal[
            "logreg", "absolut", "svm", "mlp", "nn", "randforest", None, "sum", "mc"
        ],
        model: Path = None,
    ) -> list:
        """
        Return:
            TODO return the same type
        """
        if search_type == "sum" or search_type == "mc" or search_type == None:
            ################ TODO REMOVE AFTER BENCHMARKING ################
            # load model
            if os.path.exists(model):
                printl(f"Loading model from {model}")
                with open(model, "rb") as f:
                    smodel = pickle.load(f)
            else:
                raise ValueError(f"{model} does not exist.")
            # predict
            proba = smodel.predict_proba(F)

            class_names = smodel.classes_

            output = [{cls: p for cls, p in zip(class_names, row)} for row in proba]
            return output
            ################ TODO REMOVE AFTER BENCHMARKING ################
        else:
            if search_type == "logreg":
                return self.fes_logistic_regression_model.predict(F)
            elif search_type == "absolut":
                return self.fes_absolut_search_engine.predict(F)
            elif search_type == "svm":
                return self.fes_svm_model.predict(F)
            elif search_type == "mlp":
                return self.fes_pocket_mlp_classifier_model.predict(F)
            elif search_type == "randforest":
                return self.fes_random_forest_model.predict(F)
            # elif search_type == "nn":
            #     return self._nn_predict(F)
            else:
                printl("Use valid earch option. [logreg, absolut]")
                return []

    def _load_databasefile(self, db_path: Path):
        """ """
        return pd.read_csv(db_path, sep="\t")



    

    

    def _counter_to_formula(self, counter: Counter) -> str:
        return "".join(f"{count}{key}" for key, count in sorted(counter.items()))



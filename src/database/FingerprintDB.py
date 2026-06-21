import os
from pathlib import Path
from src.utils.geometric.point_distance import point_distance
from typing import Literal
import pickle
import numpy as np
import pandas as pd
from src.io.printl import printl
from src.utils.constants import (
    FINGERPRINT_DB,
    MODEL,
    AA_ORDER,
    FINGERPRINT_RADIUS,
    MIN_DIST_ARTIFICAL_SAMPLES,
    AMOUNT_ARTIFICAL_SAMPLES
)
from sklearn.linear_model import LogisticRegression
from collections import Counter
from src.utils.smiles import SMILES
from src.filter.apply_blacklist import apply_blacklist_build
from src.filter.apply_whitelist import apply_whitelist_build
from src.utils.geometric.calculate_geometric_centers import calculate_geometric_centers
from src.fingerprint.create_fingerprint import create_fingerprint
import uuid
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from src.models.logistic_regression_model import LogisticRegressionModel
from src.models.absolut_search import AbsolutSearchEngine
from src.models.svm import SupportVectorMachine
from src.models.mlp_classifier import CustomMLPClassifier
from src.models.random_forest import CustomRandomForestClassifier

import warnings
warnings.simplefilter("ignore", PDBConstructionWarning)
warnings.filterwarnings(
    "ignore",
    message="X does not have valid feature names"
)

class FingerprintDB:
    def __init__(self, db_directory : Path = None):
        """
        """
        
        self.db_path = Path(os.path.join(str(db_directory), "fingerprint.tsv"))

        # models
        self.logistic_regression_model_path = os.path.join(str(db_directory), "logistic_regression_model.pkl")
        self.mlp_classifier_model_path = os.path.join(str(db_directory), "mlp_classifier_model.pkl")
        self.svm_model_path = os.path.join(str(db_directory), "svm_model.pkl")
        self.random_forest_path = os.path.join(str(db_directory), "randomforest_model.pkl")

        # load db
        if os.path.exists(self.db_path):
            printl(f"Loading fingerprint database from {self.db_path} ...")
            self.db = self._load_databasefile(self.db_path)

        # load models TODO no reason to load all models at once

        # load absolut search engine
        if os.path.exists(self.db_path):
            self.absolut_search_engine = AbsolutSearchEngine().load_file(self.db_path)

        
        # logistic regr model
        if os.path.exists(self.logistic_regression_model_path):
            self.logistic_regression_model = LogisticRegressionModel().load(self.logistic_regression_model_path)

        # mlp classifier model
        if os.path.exists(self.mlp_classifier_model_path):
            self.mlp_classifier_model = CustomMLPClassifier().load(self.mlp_classifier_model_path)

        # svm model
        if os.path.exists(self.svm_model_path):
            self.svm_model = SupportVectorMachine().load(self.svm_model_path)

        # random forest
        if os.path.exists(self.random_forest_path):
            printl(f"Loading model from {self.random_forest_path}")
            with open(self.random_forest_path, "rb") as f:
                self.nn_model = pickle.load(f)


    def load(self, db_path : Path = FINGERPRINT_DB, model_path : Path = MODEL):
        """
        TODO unused till here
        """
        if db_path != None:
            self.db_path = db_path
            if os.path.exists(self.db_path):
                printl(f"Loading fingerprint database from {self.db_path} ...")
                self.db = self._load_databasefile(self.db_path)

        if model_path != None:
            self.model_path = model_path
            if os.path.exists(self.model_path):
                printl(f"Loading model from {self.model_path}")
                with open(self.model_path, "rb") as f:
                    self.model = pickle.load(f)   

    def save(self, save_models : bool = True):
        """
        """
        # save db
        self.db.to_csv(self.db_path, sep="\t", index=None)

        if save_models:
            # save models
            # logreg model
            self.logistic_regression_model.save(self.logistic_regression_model_path)
            
            # svm model
            self.svm_model.save(self.svm_model_path)

            # mlp
            self.mlp_classifier_model.save(self.mlp_classifier_model_path)
                
            # random forest
            self.random_forest_model.save(self.random_forest_path)
        
    def save_fingerprint_tsv(self, path : Path):
        """"""
        self.db.to_csv(path, sep="\t", index = None)

    def search(self, F : dict, search_type : Literal["logreg", "absolut", "svm", "mlp", "nn", "randforest", None, "sum", "mc"], model : Path = None) -> list:
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

            output = [
                {cls: p for cls, p in zip(class_names, row)}
                for row in proba
            ]
            return output
            ################ TODO REMOVE AFTER BENCHMARKING ################
        else:
            if search_type == "logreg":
                return self.logistic_regression_model.predict(F)
            elif search_type == "absolut":
                return self.absolut_search_engine.predict(F)
            elif search_type == "svm":
                return self.svm_model.predict(F)
            elif search_type == "mlp":
                return self.mlp_classifier_model.predict(F)
            elif search_type == "randforest":
                return self.random_forest_model.predict(F)
            # elif search_type == "nn":
            #     return self._nn_predict(F)
            else:
                printl("Use valid earch option. [logreg, absolut]")
                return []
        
    def build_db(self, input_directory : Path, f_radius : float = FINGERPRINT_RADIUS, train_models : bool = True, extend_background_samples : bool = False):
        """
        """
        # build db table
        self.db = self._build_db_table(input_directory, f_radius, extend_background_samples)
        self.absolut_search_engine = AbsolutSearchEngine().load_dataframe(self.db_path)


        if train_models:

            X = self.db.drop(columns=['formula', 'id', 'smiles', 'res_name'])
            X_np = X.to_numpy()

            Y = self.db["formula"]
            # train_log_reg
            self.logistic_regression_model = LogisticRegressionModel()
            self.logistic_regression_model.train(X_np, Y)

            # train svm
            self.svm_model = SupportVectorMachine()
            self.svm_model.train(X_np, Y)

            # train mlp classifier
            self.mlp_classifier_model = CustomMLPClassifier()
            self.mlp_classifier_model.train(X_np, Y)

            # train random forest model
            self.random_forest_model = CustomRandomForestClassifier()
            self.random_forest_model.train(X_np, Y)

    def _load_databasefile(self, db_path : Path):
        """
        """
        return pd.read_csv(db_path, sep="\t")

    def _extract_hetatm_residues(self, pdb_file, exclude_water=True):
        """
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("struct", pdb_file)

        results = {}
        for model in structure:
            for chain in model:
                for residue in chain:

                    hetflag, resseq, icode = residue.id

                    if not hetflag.startswith("H_"):
                        continue

                    if exclude_water and residue.resname in {"HOH", "WAT"}:
                        continue

                    atoms = []
                    for atom in residue:
                        element = atom.element.strip()
                        x, y, z = atom.coord
                        atoms.append((element, float(x), float(y), float(z)))

                    unique_identifier = uuid.uuid4()

                    results[unique_identifier] = {
                        "res_name": residue.resname,
                        "res_id": resseq,
                        "chain": chain.id,
                        "atoms": atoms
                    }

        return dict(results)
    
    def _build_db_table(self, input_directory : Path, f_radius : float, extend_background_samples : bool = False):
        """
        """
        db_df = pd.DataFrame({
                "id" : [],
                "res_name": [],
                "smiles": [],
                "formula": [],
                "ALA": [],
                "ARG": [],
                "ASN": [],
                "ASP": [],
                "CYS": [],
                "GLN": [],
                "GLU": [],
                "GLY": [],
                "HIS": [],
                "ILE": [],
                "LEU": [],
                "LYS": [],
                "MET": [],
                "PHE": [],
                "PRO": [],
                "SER": [],
                "THR": [],
                "TRP": [],
                "TYR": [],
                "VAL": []
            })
        
        counter87 = 0

        all = len(list(os.listdir(input_directory)))

        # Iterate over each structure in input_structure_dir
        for structure in os.listdir(input_directory):
            printl(f"{counter87}/{all}")
            counter87 += 1
            
            structure_path = os.path.join(input_directory, structure)
            if not os.path.exists(structure_path):
                if True: #TODO change to verbose
                    printl(f"{structure_path} does not exist.")
                continue

            # load structure and identify all cofactors exept the ones from blacklist
            hetatms = self._extract_hetatm_residues(structure_path)
            printl(f"Structure {structure} has {len(hetatms)} cofactors before filtering.")

            cofactors, blacklisted = apply_blacklist_build(hetatms)
            printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering with blacklist.")

            cofactors, _ = apply_whitelist_build(cofactors)
            printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering.")
            # identify geometric center
            cofactors = calculate_geometric_centers(cofactors, "atoms")
            #print(cofactors)
            printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering and calculating geometric centers.")
            
            blacklisted = calculate_geometric_centers(blacklisted, "atoms")            

            # for each cofactor left
            for c in cofactors:
                # create fingerprint
                F = create_fingerprint(structure_path, cofactors[c]["geometric_center"], f_radius)
                
                F_dict = {i: F[AA_ORDER.index(i)] for i in AA_ORDER}

                # create formula
                atoms = Counter(atom[0].upper() for atom in cofactors[c]["atoms"]
                    if atom[0].upper() in {"FE", "S"})
                
                atoms_lst = [atom[0] for atom in cofactors[c]["atoms"]]

                if "NI" in atoms_lst or "N" in atoms_lst:
                    formula  = "active_site"
                else:                
                    formula = self._counter_to_formula(atoms)

                # skip anything thats not 3/4FE3/4S
                if formula not in ("3FE4S", "4FE3S", "4FE4S", "active_site"): # TODO open at some point for other fes clusters
                    continue

                F_dict["formula"] = formula
                
                try:
                    F_dict["smiles"] = [SMILES[formula]] #TODO change to smiles
                except:
                    F_dict["smiles"] = "UwU"

                F_dict["id"] = [c]
                F_dict["res_name"] = [cofactors[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)

            # test none class
            for c in blacklisted:
                # create Fingerprint
                F = create_fingerprint(structure_path, blacklisted[c]["geometric_center"])
                F_dict = {i: F[AA_ORDER.index(i)] for i in AA_ORDER}
                F_dict["formula"] = "protein"
                F_dict["smiles"] = "None"
                F_dict["id"] = [c]
                F_dict["res_name"] = [blacklisted[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)

            # if extend_protein_samples
            if extend_background_samples:
                """
                Strategy B:
                1. get boundaries of protein
                2. create distributions from boundaries
                3. generate points using these distributions (always checking if they are to close to a cofactor geometric center)
                4. calculate fingerprint and add
                """
                # 1.
                x_min, x_max, y_min, y_max, z_min, z_max = self._get_protein_boundaries(structure_path)

                mins = np.array([x_min, y_min, z_min])
                maxs = np.array([x_max, y_max, z_max])

                # 2.
                center = (mins + maxs) / 2
                sigma = (maxs - mins) / 6

                # 3.
                samples = np.random.normal(loc=center, scale=sigma, size=(AMOUNT_ARTIFICAL_SAMPLES, 3))
                mask = np.all((samples >= mins) & (samples <= maxs), axis=1)
                samples = samples[mask]

                counter89 = 0

                for sample in samples:
                    # check distance to cofactor geometric centers
                    if all(point_distance(sample[0], sample[1], sample[2], cofactors[o]["geometric_center"][0], cofactors[o]["geometric_center"][1], cofactors[o]["geometric_center"][2]) >= MIN_DIST_ARTIFICAL_SAMPLES for o in cofactors.keys()):
                        F = create_fingerprint(structure_path, sample)
                        F_dict = {i: F[AA_ORDER.index(i)] for i in AA_ORDER}
                        F_dict["formula"] = "protein"
                        F_dict["smiles"] = "None"
                        F_dict["id"] = ["None"]
                        F_dict["res_name"] = ["artificial"]

                        db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)
                        counter89 += 1
                printl(f"Added {counter89} artificial samples.")
                


            printl(f"Structrue {structure} complete ...")

        return db_df

    def _counter_to_formula(self, counter: Counter) -> str:
        return "".join(f"{count}{key}" for key, count in sorted(counter.items()))
    
    def _get_protein_boundaries(self, structure_path : Path) -> tuple:
        """
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('s', structure_path)

        coords = np.array([atom.coord for atom in structure.get_atoms()])
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)

        x_min, x_max, y_min, y_max, z_min, z_max = mins[0], maxs[0], mins[1], maxs[1], mins[2], maxs[2]
        return x_min, x_max, y_min, y_max, z_min, z_max
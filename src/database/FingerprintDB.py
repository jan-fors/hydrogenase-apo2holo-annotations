import os
from pathlib import Path
from typing import Literal
import pickle
import pandas as pd
from src.io.printl import printl
from src.utils.constants import (
    FINGERPRINT_DB,
    MODEL,
    MAX_ITER,
    SOLVER
)
from sklearn.linear_model import LogisticRegression
from collections import Counter
from src.utils.smiles import SMILES
from src.filter.apply_blacklist import apply_blacklist_build
from src.filter.apply_whitelist import apply_whitelist_build
from src.utils.calculate_geometric_centers import calculate_geometric_centers
from src.fingerprint.create_fingerprint import create_fingerprint
import uuid
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.simplefilter("ignore", PDBConstructionWarning)
import warnings

warnings.filterwarnings(
    "ignore",
    message="X does not have valid feature names"
)

class FingerprintDB:
    def __init__(self, db_directory : Path = None):
        """
        """
        
        self.db_path = os.path.join(str(db_directory), "fingerprint.tsv")
        self.model_path = os.path.join(str(db_directory), "model.pkl")

        if os.path.exists(self.db_path):
            printl(f"Loading fingerprint database from {self.db_path} ...")
            self.db = self._load_databasefile(self.db_path)

        if os.path.exists(self.model_path):
            printl(f"Loading model from {self.model_path}")
            with open(self.model_path, "rb") as f:
                self.model = pickle.load(f)     
    
    def load(self, db_path : Path = FINGERPRINT_DB, model_path : Path = MODEL):
        """
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

    def save(self, db_path : Path = FINGERPRINT_DB, model_path : Path = MODEL):
        """
        """
        if db_path != None:
            self.db.to_csv(db_path, sep="\t", index=None)

        if model_path != None:
            with open(model_path, "wb") as f:
                pickle.dump(self.model, f)
    
    def search(self, F : dict, search_type : Literal["logreg", "absolut"]) -> list:
        """
        """
        if search_type == "logreg":
            return self._predict(F)
        elif search_type == "absolut":
            return self._absolut_search(F)
        else:
            printl("Use valid earch option. [logreg, absolut]")
            return []

    def get_model(self):
        return self.model

    def build_db(self, input_directory : Path):
        """
        """
        # build db table
        self.db = self._build_db_table(input_directory)

        # train_log_reg
        self.model = self._train_log_reg(self.db)

    def _train_log_reg(self, database : pd.DataFrame):
        """
        
        """
        # database = database[~database['formula'].str.contains('NI', na=False)]
        # database = database[~database['formula'].str.contains('N', na=False)]
        # database = database[~database['formula'].str.contains('O', na=False)]

        X = database.drop(columns=['formula', 'id', 'smiles', 'res_name'])
        y = database["formula"]

        model = LogisticRegression(solver=SOLVER, max_iter=MAX_ITER, class_weight="balanced")
        model.fit(X, y)

        return model

    def _predict(self, F : dict) -> list:
        """
        """
        #X = [list(F.values())]
        X = [[F[name] for name in self.model.feature_names_in_]]
        proba = self.model.predict_proba(X)
        # proba_max = float(proba[0].max())
        # if proba_max < 0.6:
        #     return None
        return [proba]

    def _absolut_search(self, F : dict) -> list:
        """
        """
        subset = self.db[self.db["ALA"] == F["ALA"]]
        subset = subset[subset["ARG"] == F["ARG"]]
        subset = subset[subset["ASN"] == F["ASN"]]
        subset = subset[subset["ASP"] == F["ASP"]]
        subset = subset[subset["CYS"] == F["CYS"]]
        subset = subset[subset["GLN"] == F["GLN"]]
        subset = subset[subset["GLU"] == F["GLU"]]
        subset = subset[subset["GLY"] == F["GLY"]]
        subset = subset[subset["HIS"] == F["HIS"]]
        subset = subset[subset["ILE"] == F["ILE"]]
        subset = subset[subset["LEU"] == F["LEU"]]
        subset = subset[subset["LYS"] == F["LYS"]]
        subset = subset[subset["MET"] == F["MET"]]
        subset = subset[subset["PHE"] == F["PHE"]]
        subset = subset[subset["PRO"] == F["PRO"]]
        subset = subset[subset["SER"] == F["SER"]]
        subset = subset[subset["THR"] == F["THR"]]
        subset = subset[subset["TRP"] == F["TRP"]]
        subset = subset[subset["TYR"] == F["TYR"]]
        subset = subset[subset["VAL"] == F["VAL"]]

        return subset["formula"].to_list()

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
    
    def _build_db_table(self, input_directory : Path):
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
        
        # Iterate over each structure in input_structure_dir
        for structure in os.listdir(input_directory):
            
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
            printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering and calculating geometric centers.")
            
            blacklisted = calculate_geometric_centers(blacklisted, "atoms")            

            # for each cofactor left
            for c in cofactors:
                # create fingerprint
                F = create_fingerprint(structure_path, cofactors[c]["geometric_center"])
                
                # create formula
                atoms = Counter(atom[0].upper() for atom in cofactors[c]["atoms"]
                    if atom[0].upper() in {"FE", "S"})
                
                atoms_lst = [atom[0] for atom in cofactors[c]["atoms"]]

                if "NI" in atoms_lst or "N" in atoms_lst:
                    formula  = "active_site"
                else:                
                    formula = self._counter_to_formula(atoms)
               
                F["formula"] = formula
                
                try:
                    F["smiles"] = [SMILES[formula]] #TODO change to smiles
                except:
                    F["smiles"] = "UwU"

                F["id"] = [c]
                F["res_name"] = [cofactors[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F)], ignore_index=True)

            # test none class
            for c in blacklisted:
                # create Fingerprint
                F = create_fingerprint(structure_path, blacklisted[c]["geometric_center"])

                F["formula"] = "protein"
                F["smiles"] = "OwO"
                F["id"] = [c]
                F["res_name"] = [blacklisted[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F)], ignore_index=True)


            printl(f"Structrue {structure} complete ...")

        return db_df

    def _counter_to_formula(self, counter: Counter) -> str:
        return "".join(f"{count}{key}" for key, count in sorted(counter.items()))
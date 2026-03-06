from pathlib import Path
import pandas as pd
import os
import uuid
import argparse
from src.io.printl import printl
import warnings
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter("ignore", PDBConstructionWarning)
from collections import defaultdict
from src.filter.apply_blacklist import apply_blacklist
from src.filter.apply_whitelist import apply_whitelist
from src.utils.calculate_geometric_centers import calculate_geometric_centers
from src.fingerprint.create_fingerprint import create_fingerprint

def build_fingerprint_db(input_structure_dir : Path, output_dir : Path):
    """
    TODO might change this later to sqlite db but at this moment dataframe is okay.
    """
    db_df = pd.DataFrame({
        "id" : [],
        "res_name": [],
        "smiles_tmp": [],
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
    for structure in os.listdir(input_structure_dir):
        
        structure_path = os.path.join(input_structure_dir, structure)
        if not os.path.exists(structure_path):
            if True: #TODO change to verbose
                printl(f"{structure_path} does not exist.")
            continue

        # load structure and identify all cofactors exept the ones from blacklist
        hetatms = _extract_hetatm_residues(structure_path)

        cofactors = apply_blacklist(hetatms)

        cofactors = apply_whitelist(cofactors)

        # identify geometric center
        cofactors = calculate_geometric_centers(cofactors, "atoms")
        
        # for each cofactor left
        for c in cofactors:
            # create fingerprint
            F = create_fingerprint(structure_path, cofactors[c]["geometric_center"])
            
            # append to database
            smiles_tmp = ""
            for atom in cofactors[c]["atoms"]:
                smiles_tmp += "-"+atom[0]

            F["smiles_tmp"] = [smiles_tmp] #TODO change to smiles
            F["id"] = [c]
            F["res_name"] = [cofactors[c]["res_name"]]

            db_df = pd.concat([db_df, pd.DataFrame(F)], ignore_index=True)

        print(f"Structrue {structure} complete ...")

    os.makedirs(output_dir, exist_ok=True)
    db_df.to_csv(os.path.join(output_dir, "fingerprint_db.tsv"),sep="\t")

def _extract_hetatm_residues(pdb_file, exclude_water=True):
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Path to the structure folder.")
    parser.add_argument("output_dir", help="Path to the structure files.")
    args = parser.parse_args()
    build_fingerprint_db(args.input_dir, args.output_dir)
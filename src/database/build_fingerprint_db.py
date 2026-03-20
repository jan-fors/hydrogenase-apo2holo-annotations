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
from src.filter.apply_blacklist import apply_blacklist_build
from src.filter.apply_whitelist import apply_whitelist_build
from src.utils.calculate_geometric_centers import calculate_geometric_centers
from src.fingerprint.create_fingerprint import create_fingerprint

def build_fingerprint_db(input_structure_dir : Path, output_dir : Path, sample : bool = False):
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
    counter = 0
    # Iterate over each structure in input_structure_dir
    for structure in os.listdir(input_structure_dir):
        
        structure_path = os.path.join(input_structure_dir, structure)
        if not os.path.exists(structure_path):
            if True: #TODO change to verbose
                printl(f"{structure_path} does not exist.")
            continue

        # load structure and identify all cofactors exept the ones from blacklist
        hetatms = _extract_hetatm_residues(structure_path)
        printl(f"Structure {structure} has {len(hetatms)} cofactors before filtering.")

        cofactors = apply_blacklist_build(hetatms)
        printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering with blacklist.")

        cofactors = apply_whitelist_build(cofactors)
        printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering.")
        # identify geometric center
        cofactors = calculate_geometric_centers(cofactors, "atoms")
        printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering and calculating geometric centers.")
        
        # for each cofactor left
        for c in cofactors:
            # create fingerprint
            F = create_fingerprint(structure_path, cofactors[c]["geometric_center"])
            
            # append to database
            smiles_tmp = ""
            for atom in sorted(cofactors[c]["atoms"], key=lambda x: x[0]):
                smiles_tmp += "-"+atom[0]

            F["smiles_tmp"] = [smiles_tmp] #TODO change to smiles
            F["id"] = [c]
            F["res_name"] = [cofactors[c]["res_name"]]

            print(F)
            db_df = pd.concat([db_df, pd.DataFrame(F)], ignore_index=True)


        print(f"Structrue {structure} complete ...")
        counter += 1
        if sample and counter >= 5:
            break

    # create row that combines all aminoacid rows into one string for deduplication
    #db_df["amino_acids"] = db_df[["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]].apply(lambda row: "".join([f"{col}:{row[col]}" for col in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]]), axis=1)
    # deduplicate database by all columns except id and res_name and smiles_tmp
    #db_df = db_df.drop_duplicates(subset=["amino_acids"])
    #db_df = db_df.drop(columns=["amino_acids"])
    #print(f"Database size after deduplication: {db_df.shape[0]}")

    os.makedirs(output_dir, exist_ok=True)
    db_df.to_csv(os.path.join(output_dir, "fingerprint_db.tsv"),sep="\t", index=None)

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
    parser.add_argument("--sample", action="store_true", help="Whether to sample the database for testing.")
    args = parser.parse_args()
    build_fingerprint_db(args.input_dir, args.output_dir, args.sample)
from pathlib import Path
from typing import Literal
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from apo2holo.io.writers.printl import printl
from apo2holo.utils.geometric.point_distance import point_distance
from apo2holo.utils.geometric.calculate_geometric_centers import calculate_geometric_centers
from apo2holo.utils.geometric.generate_sphere_points import random_points_in_sphere, fibonacci_sphere
from apo2holo.filter.apply_blacklist import apply_blacklist_build
from apo2holo.filter.apply_whitelist import apply_whitelist_build
from apo2holo.fingerprint.create_fingerprint import create_fingerprint
from apo2holo.utils.formula.smiles import SMILES
from apo2holo.pdb.pdb_handler import extract_hetatm_residues, get_protein_boundaries
import numpy as np
from apo2holo.utils.formula.counter_to_formula import counter_to_formula

import traceback

def build_fingerprint_tsv(input_directory: Path,
        f_radius: float,
        extend_background_samples: bool = False,
        cofactor_augmentation: bool = False,
        threads: int = 1,
        fingerprint_types : list = [],
        min_dist_art_samples : float = 0.0,
        n_art_samples : int = 0,
        strategy : Literal["fibonacci", "random"] = "random",
        aug_radius : float = 0.0,
        n_augs : int = 0) -> pd.DataFrame:
    """
    """
    structures = os.listdir(input_directory)
    
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                _process_structure,
                file,
                input_directory,
                f_radius,
                fingerprint_types,
                extend_background_samples,
                min_dist_art_samples,
                n_art_samples,
                strategy,
                cofactor_augmentation,
                aug_radius,
                n_augs
            ): file
            for file in structures
        }
        for future in as_completed(futures):
            file = futures[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
            except Exception as e:
                printl(f"Failed on {file}: {e}")

    final_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    return final_df

def build_and_safe_fingerprint_tsv(
        input_directory: Path,
        out_file: Path,
        f_radius: float,
        extend_background_samples: bool = False,
        cofactor_augmentation: bool = False,
        threads: int = 1,
        fingerprint_types : list = [],
        min_dist_art_samples : float = 0.0,
        n_art_samples : int = 0,
        strategy : Literal["fibonacci", "random"] = "random",
        aug_radius : float = 0.0,
        n_augs : int = 0):
    """
    """
    structures = os.listdir(input_directory)
        
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(
                _process_structure,
                file,
                input_directory,
                f_radius,
                fingerprint_types,
                extend_background_samples,
                min_dist_art_samples,
                n_art_samples,
                strategy,
                cofactor_augmentation,
                aug_radius,
                n_augs
            ): file
            for file in structures
        }
        for future in as_completed(futures):
            file = futures[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
            except Exception as e:
                printl(f"Failed on {file}: {e}")
                traceback.print_exc()
                exit(0)

    final_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    final_df.to_csv(out_file, sep="\t", index=None)

def _process_structure(
        structure: str,
        input_directory: Path,
        f_radius,
        fingerprint_types : list = [],
        extend_background_samples: bool = False,
        min_dist_art_samples: float = 0.0,
        n_art_samples : int = 0,
        strategy : Literal["fibonacci", "random"] = "random",
        cofactor_augmentation: bool = False,
        aug_radius : float = 0.0,
        n_augs : int = 0
    ) -> pd.DataFrame:
        """ """
        db_df = pd.DataFrame()
        structure_path = os.path.join(input_directory, structure)
        if not os.path.exists(structure_path):
            printl(f"{structure_path} does not exist.")

        # load structure and identify all cofactors exept the ones from blacklist
        hetatms = extract_hetatm_residues(structure_path)
        
        printl(f"Structure {structure} has {len(hetatms)} cofactors before filtering.")

        cofactors, blacklisted = apply_blacklist_build(hetatms)

        printl(
            f"Structure {structure} has {len(cofactors)} cofactors after filtering with blacklist."
        )

        cofactors, _ = apply_whitelist_build(cofactors)

        printl(f"Structure {structure} has {len(cofactors)} cofactors after filtering.")

        # identify geometric center
        cofactors = calculate_geometric_centers(cofactors, "atoms")

        # print(cofactors)
        printl(
            f"Structure {structure} has {len(cofactors)} cofactors after filtering and calculating geometric centers."
        )

        blacklisted = calculate_geometric_centers(blacklisted, "atoms")

        # for each cofactor left
        for c in cofactors:
            if cofactor_augmentation:
                g_center = cofactors[c]["geometric_center"]
                a_centers = random_points_in_sphere(g_center, aug_radius, n_augs)
                centers = [(g_center, 0.0)]
                for i in a_centers:
                    centers.append((i, point_distance(g_center[0], g_center[1], g_center[2], i[0], i[1], i[2])))

            else:
                centers = [(cofactors[c]["geometric_center"], 0.0)]

            for center in centers:
                # create fingerprint
                F = create_fingerprint(structure_path, center[0], f_radius, fingerprint_types)

                F_dict = {i: value for i, value in enumerate(F)}
                F_dict["aug_dist"] = center[1]
                F_dict["structure"] = structure
                F_dict["type"] = "cofactor"

                # create formula
                atoms = Counter(
                    atom[0].upper()
                    for atom in cofactors[c]["atoms"]
                    if atom[0].upper() in {"FE", "S"}
                )

                atoms_lst = [atom[0] for atom in cofactors[c]["atoms"]]

                if "NI" in atoms_lst or "N" in atoms_lst:
                    formula = "active_site"
                else:
                    formula = counter_to_formula(atoms)

                # skip anything thats not 3/4FE3/4S
                if formula not in (
                    "3FE4S",
                    "4FE3S",
                    "4FE4S",
                    "active_site",
                ):  # TODO open at some point for other fes clusters
                    continue

                F_dict["formula"] = formula

                try:
                    F_dict["smiles"] = [SMILES[formula]]  # TODO change to smiles
                except:
                    F_dict["smiles"] = "UwU"

                F_dict["id"] = [c]
                F_dict["res_name"] = [cofactors[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)

        # test none class
        for c in blacklisted:
            if cofactor_augmentation:
                g_center = blacklisted[c]["geometric_center"]
                a_centers = random_points_in_sphere(g_center, aug_radius, n_augs)

                centers = [(g_center, 0.0)]

                for i in a_centers:
                    centers.append((i, point_distance(g_center[0], g_center[1], g_center[2], i[0], i[1], i[2])))
            else:
                centers = [(cofactors[c]["geometric_center"], 0.0)]

            for center in centers:

                # create Fingerprint
                F = create_fingerprint(structure_path, center[0], f_radius, fingerprint_types)
                F_dict = {i: value for i, value in enumerate(F)}
                F_dict["aug_dist"] = center[1]
                F_dict["structure"] = structure
                F_dict["type"] = "blacklisted cofactor"
                F_dict["formula"] = "protein"
                F_dict["smiles"] = "None"
                F_dict["id"] = [c]
                F_dict["res_name"] = [blacklisted[c]["res_name"]]

                db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)

        # if extend_protein_samples TODO
        if extend_background_samples:
            """
            Strategy B:
            1. get boundaries of protein
            2. create distributions from boundaries
            3. generate points using these distributions (always checking if they are to close to a cofactor geometric center)
            4. calculate fingerprint and add
            """
            # opt A
            if strategy == "random":
                # 1.
                x_min, x_max, y_min, y_max, z_min, z_max = get_protein_boundaries(
                    structure_path
                )

                mins = np.array([x_min, y_min, z_min])
                maxs = np.array([x_max, y_max, z_max])

                # 2.
                center = (mins + maxs) / 2
                sigma = (maxs - mins) / 6

                # 3.
                samples = np.random.normal(
                    loc=center, scale=sigma, size=(n_art_samples, 3)
                )
                mask = np.all((samples >= mins) & (samples <= maxs), axis=1)
                samples = samples[mask]

                counter89 = 0

                for sample in samples:
                    # check distance to cofactor geometric centers
                    if all(
                        point_distance(
                            sample[0],
                            sample[1],
                            sample[2],
                            cofactors[o]["geometric_center"][0],
                            cofactors[o]["geometric_center"][1],
                            cofactors[o]["geometric_center"][2],
                        )
                        >= min_dist_art_samples
                        for o in cofactors.keys()
                    ):
                        F = create_fingerprint(structure_path, center[0], f_radius, fingerprint_types)
                        F_dict = {i: value for i, value in enumerate(F)}
                        F_dict["aug_dist"] = 0.0
                        F_dict["structure"] = structure
                        F_dict["type"] = "random background point"
                        F_dict["formula"] = "protein"
                        F_dict["smiles"] = "None"
                        F_dict["id"] = ["None"]
                        F_dict["res_name"] = ["artificial"]

                        db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)
                        counter89 += 1
                printl(f"Added {counter89} artificial samples.")
            
            else:
                for c in cofactors:
                    g_center = cofactors[c]["geometric_center"]

                    background_centers = fibonacci_sphere(g_center, f_radius, n_art_samples)
                    
                    for bc in background_centers:
                        F = create_fingerprint(structure_path, center[0], f_radius, fingerprint_types)
                        F_dict = {i: value for i, value in enumerate(F)}
                        F_dict["aug_dist"] = 0.0
                        F_dict["structure"] = structure
                        F_dict["type"] = "fibonacci sphere cofactor"
                        F_dict["formula"] = "protein"
                        F_dict["smiles"] = "None"
                        F_dict["id"] = ["None"]
                        F_dict["res_name"] = ["artificial"]

                        db_df = pd.concat([db_df, pd.DataFrame(F_dict)], ignore_index=True)

        printl(f"Structrue {structure} complete ...")
        return db_df
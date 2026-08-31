import os
from pathlib import Path
from typing import Literal
from src.database.structure.handler import identify_possible_active_site_pockets
import numpy as np
from src.utils.geometric.calculate_geometric_centers import (
    calculate_center,
)
from src.fingerprint.create_fingerprint import create_fingerprint
from src.io.writers.printl import printl, print_probabilities
from src.utils.geometric.point_distance import point_distance
from collections import Counter
from src.database.FingerprintDB import FingerprintDB
from typing import List


######################## HELPER FUNCTIONS ########################


def _counter_to_probs(counter):
    total = sum(counter.values())
    if total == 0:
        return {k: 0 for k in counter}
    return {k: v / total for k, v in counter.items()}










def _identify_fes_cluster(best_hits: dict) -> List[int]:
    """
    best hits contains:
    key -> cluster number
    item -> name of cofactor
    """
    keep = []
    for h in best_hits:
        if best_hits[h] != "protein":
            keep.append(h)

    return keep


def _sort_into_hyd_format(
    active_site_mass_center: tuple, fes_cluster_mass_center: dict
) -> tuple:
    """ """
    proximal_key = None
    medial_key = None
    distal_key = None

    distances = []

    for k in fes_cluster_mass_center.keys():
        dist = point_distance(
            active_site_mass_center[0],
            active_site_mass_center[1],
            active_site_mass_center[2],
            fes_cluster_mass_center[k][0],
            fes_cluster_mass_center[k][1],
            fes_cluster_mass_center[k][2],
        )
        distances.append((k, dist))

    # sort distances by dist
    sorted_distances = sorted(distances, key=lambda x: x[1])

    try:
        proximal_key = sorted_distances[0][0]
    except IndexError:
        printl("Not enough iron-sulfur clusters identified for proximal annotation")

    try:
        medial_key = sorted_distances[1][0]
    except IndexError:
        printl("Not enough iron-sulfur clusters identified for medial annotation")

    try:
        distal_key = sorted_distances[2][0]
    except IndexError:
        printl("Not enough iron-sulfur clusters identified for distal annotation")

    return proximal_key, medial_key, distal_key


######################## HELPER FUNCTIONS ########################
def fit_into_hyd_structure(
    active_site_mass_center: tuple, distances: List, best_hits: dict, pockets: dict
) -> tuple:
    """ """

    fes_cluster_keys = _identify_fes_cluster(best_hits)
    fes_cluster_mass_centers = {}
    for k in fes_cluster_keys:
        fes_cluster_mass_centers[k] = calculate_center(np.array(pockets[k]))

    proximal_key, medial_key, distal_key = _sort_into_hyd_format(
        active_site_mass_center, fes_cluster_mass_centers
    )

    return proximal_key, medial_key, distal_key














def _predict_fes_prob_by_pockets(
    fingerprint_db_path: Path,
    pockets: dict,
    search_type: str,
    input_path: Path,
    model: Path = None,
    f_radius: float = None,
):
    # Initialize fingerprintDB
    printl("Initializing fes fingerprint database...")
    fes_fingerprint_db_path = fingerprint_db_path / Path(
        "fingerprintDB"
    )  # TODO is this needed?
    fes_fingerprintDB = FingerprintDB(fes_fingerprint_db_path)

    # search each pocket against
    printl("Searching each pocket against the fes type fingerprint database ...")

    pocket_results = {}

    for key in pockets.keys():
        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(
                input_path=input_path, coordinates=pockets[key], f_radius=f_radius
            )

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = fes_fingerprintDB.fes_type_search(
                F=fingerprints, search_type=search_type, model=model
            )

            # count
            pocket_results[key] = Counter()
            pocket_results[key].update(hits)

            pocket_results[key] = _counter_to_probs(pocket_results[key])

        elif (
            search_type == "logreg_sum"
            or search_type == "mlp_sum"
            or search_type == "svm_sum"
            or search_type == "randforest_sum"
            or search_type == "sum"
        ):
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(
                input_path=input_path, coordinates=pockets[key], f_radius=f_radius
            )

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fes_fingerprintDB.fes_type_search(
                F=fingerprints, search_type=search_type_command, model=model
            )

            # count
            pocket_results[key] = _aggregate_probs(hits)

        elif (
            search_type == "logreg_mc"
            or search_type == "mlp_mc"
            or search_type == "svm_mc"
            or search_type == "randforest_mc"
            or search_type == "mc"
        ):
            # create center
            M = np.array(pockets[key])
            center = calculate_center(M)
            printl(f"Creating fingerprint for pocket {key}...")
            if f_radius != None:
                F = [create_fingerprint(input_path, center, f_radius)]
            else:
                F = [create_fingerprint(input_path, center)]

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fes_fingerprintDB.fes_type_search(
                F=F, search_type=search_type_command, model=model
            )

            pocket_results[key] = _counter_to_probs(Counter(hits[0]))


        print("="*30, "PROBS", "="*30)
        print_probabilities(pocket_results[key])
        print("="*67)

    return pocket_results

def _predict_fes_pockets(fingerprint_db_path: Path,
    pockets: dict,
    search_type: str,
    input_path: Path,
    model: Path = None,
    f_radius: float = None):
    """
    """
    # Initialize fingerprintDB
    printl("Initializing fes pocket fingerprint database...")
    fes_fingerprint_db_path = fingerprint_db_path / Path(
        "fingerprintDB"
    )  # TODO is this needed?
    fes_fingerprintDB = FingerprintDB(fes_fingerprint_db_path)

    # search each pocket against
    printl("Searching each pocket against the fes pocket fingerprint database ...")

    pocket_results = {}

    for key in pockets.keys():
        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(
                input_path=input_path, coordinates=pockets[key], f_radius=f_radius
            )

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = fes_fingerprintDB.fes_pocket_search(
                F=fingerprints, search_type=search_type, model=model
            )

            # count
            pocket_results[key] = Counter()
            pocket_results[key].update(hits)

            pocket_results[key] = _counter_to_probs(pocket_results[key])

        elif (
            search_type == "logreg_sum"
            or search_type == "mlp_sum"
            or search_type == "svm_sum"
            or search_type == "randforest_sum"
            or search_type == "sum"
        ):
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(
                input_path=input_path, coordinates=pockets[key], f_radius=f_radius
            )

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fes_fingerprintDB.fes_pocket_search(
                F=fingerprints, search_type=search_type_command, model=model
            )

            # count
            pocket_results[key] = _aggregate_probs(hits)

        elif (
            search_type == "logreg_mc"
            or search_type == "mlp_mc"
            or search_type == "svm_mc"
            or search_type == "randforest_mc"
            or search_type == "mc"
        ):
            # create center
            M = np.array(pockets[key])
            center = calculate_center(M)
            printl(f"Creating fingerprint for pocket {key}...")
            if f_radius != None:
                F = [create_fingerprint(input_path, center, f_radius)]
            else:
                F = [create_fingerprint(input_path, center)]

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fes_fingerprintDB.fes_pocket_search(
                F=F, search_type=search_type_command, model=model
            )

            pocket_results[key] = _counter_to_probs(Counter(hits[0]))


        print("="*30, "PROBS", "="*30)
        print_probabilities(pocket_results[key])
        print("="*67)

    return pocket_results

def _filter_fes_pockets_by_probability(fes_pockets, fes_probabilites):
    """"""
    res = {}
    for p in fes_probabilites.keys():
        if _get_best_hits_name(fes_probabilites[p]) == "fes_pocket":
            res[p] = fes_pockets[p]

    return res



def main(
    input_path: Path,
    out: Path,
    tmp: Path,
    boltz: bool,
    plot: bool,
    chain_dir: Path,
    structure_db_path: Path,
    search_type : Literal["mc", "sum"],
    active_site_model: Path = None,
    fes_type_model: Path = None,
    fes_pocket_model : Path = None,
    f_radius: float = None,
):
    """
    """
    # step 1: identify possible active site pockets
    chain_paths, active_site_pockets = identify_possible_active_site_pockets(
        structure_db_path=structure_db_path,
        input_path=input_path,
        out=out,
        chain_dir=chain_dir,
    )

    if chain_paths != None and active_site_pockets != None:
        # step 2: predict active site probability per pocket
        active_site_probabilites = _predict_active_site_prob_by_pockets(
            fingerprint_db_path=fingerprint_db_path,
            pockets=active_site_pockets,
            search_type=search_type,
            input_path=input_path,
            model=active_site_model,
            f_radius=f_radius,
        )

        # step 3: choose highest confidence active site
        active_site_pocket_key = _choose_highest_conf_active_site_pocket(
            active_site_probabilites
        )
        active_site_cluster_points = active_site_pockets[active_site_pocket_key]

        # step 4: identify possible fes pockets
        fes_pockets = _identify_possible_fes_pockets(
            structure_db_path=structure_db_path,
            chain_paths=chain_paths,
            input_path=input_path,
            out=out,
            chain_dir=chain_dir,
        )

        if fes_pockets != None:

            # step 5: check if clusters overlap with active_site cluster
            active_site_cluster_center = calculate_center(np.array(active_site_cluster_points))

            valid_keys = []
            for pocket_key in fes_pockets.keys():
                pocket_center = calculate_center(np.array(fes_pockets[pocket_key]))

                dist = point_distance(
                    active_site_cluster_center[0],
                    active_site_cluster_center[1],
                    active_site_cluster_center[2],
                    pocket_center[0],
                    pocket_center[1],
                    pocket_center[2],
                )
                if dist > NN_CLUSTERING_RADIUS:
                    valid_keys.append(pocket_key)
            for key in list(
                fes_pockets.keys()
            ):  # list() needed since you're modifying while iterating
                if key not in valid_keys:
                    del fes_pockets[key]

            # step 5.2: check whether the pockets are really fes pockets using the fes pocket model
            fes_pocket_probabilites = _predict_fes_pockets(
                fingerprint_db_path=fingerprint_db_path,
                pockets=fes_pockets,
                search_type=search_type,
                input_path=input_path,
                model=fes_pocket_model,
                f_radius=f_radius,)
            
            # step 5.3: filter fes pockets by pocket probability
            fes_pockets = _filter_fes_pockets_by_probability(fes_pockets, fes_pocket_probabilites)

            # step 6: predict each pocket
            fes_probabilities = _predict_fes_prob_by_pockets(
                fingerprint_db_path=fingerprint_db_path,
                pockets=fes_pockets,
                search_type=search_type,
                input_path=input_path,
                model=fes_type_model,
                f_radius=f_radius,
            )

            # step 7: choose highest confidence hits
            best_fes_hits_per_pocket = _choose_best_fes_per_pocket(
                prediction_per_pocket=fes_probabilities
            )

            # step 8: fit sort and fit into hyd
            cluster_dist = []
            for pocket_key in fes_pockets.keys():
                pocket_center = calculate_center(np.array(fes_pockets[pocket_key]))
                dist = point_distance(
                    active_site_cluster_center[0],
                    active_site_cluster_center[1],
                    active_site_cluster_center[2],
                    pocket_center[0],
                    pocket_center[1],
                    pocket_center[2],
                )
                cluster_dist.append((pocket_key, dist))

            sorted_list = sorted(cluster_dist, key=lambda t: t[1])

            proximal_key, medial_key, distal_key = fit_into_hyd_structure(
                active_site_mass_center=active_site_cluster_center,
                distances=sorted_list,
                best_hits=best_fes_hits_per_pocket,
                pockets=fes_pockets,
            )

            # step 9: output
            write_outputs(
                out=out,
                sorted_list=sorted_list,
                structure_path=input_path,
                as_pockets=active_site_pockets,
                as_pocket_key=active_site_pocket_key,
                as_probabilities=active_site_probabilites,
                fes_pockets=fes_pockets,
                prediction_per_pocket=fes_probabilities,
                best_hits_per_pocket=best_fes_hits_per_pocket,
                proximal_key=proximal_key,
                medial_key=medial_key,
                distal_key=distal_key,
                boltz=boltz,
                plot=plot,
            )
        else:
            printl("Ending program")
    else:
        printl("Ending program")
from pathlib import Path
from typing import List, Literal
import numpy as np

from apo2holo.io.readers.load_model import load_model
from apo2holo.io.writers.printl import printl

from apo2holo.database.structure.handler import identify_possible_active_site_pockets, identify_possible_fes_pockets
from apo2holo.models.inference.inference_handler import predict_cofactor_prob_for_pockets

from apo2holo.utils.geometric.calculate_geometric_centers import calculate_center
from apo2holo.utils.geometric.point_distance import point_distance

from apo2holo.io.writers.output_handler import write_outputs


def main(
    input_structure_path: Path,
    output_dir: Path,
    structure_db_path: Path,
    chain_dir_path: Path,
    nn_clustering_radius : float,
    fident_threshold : float,
    bits_threshold : float,
    as_model_path: Path,
    as_search_type: Literal["sum", "mc"],
    fes_pocket_model_path : Path,
    fes_pocket_search_type : Literal["sum", "mc"],
    fes_type_model_path : Path,
    fes_type_search_type : Literal["sum", "mc"],
    plot : bool,
    boltz : bool
):
    """ """
    # step 1: identify possible active site pockets
    chain_paths, active_site_pockets = identify_possible_active_site_pockets(
        structure_db_path=structure_db_path,
        input_path=input_structure_path,
        out=output_dir,
        chain_dir=chain_dir_path,
        nn_clustering_radius=nn_clustering_radius,
        fident_threshold=fident_threshold,
        bits_threshold=bits_threshold
    )

    if chain_paths == None and active_site_pockets == None:
        print("No active site pockets found.")
        exit(0)

    # step 2: predict active site probability per pocket

    # load as model
    as_model, _, fingerprint_type, f_radius = load_model(as_model_path)

    # predict active site probabilites
    active_site_probabilites = predict_cofactor_prob_for_pockets(
        structure=input_structure_path,
        model=as_model,
        fingerprint_types=fingerprint_type,
        f_radius=f_radius,
        search_type=as_search_type,
        pockets=active_site_pockets,
    )

    for key in active_site_probabilites:
        active_site_probabilites[key] = active_site_probabilites[key]["active_site"]

    # step 3: choose highest confidence active site
    active_site_pocket_key = _choose_highest_conf_active_site_pocket(
        active_site_probabilites
    )
    active_site_cluster_points = active_site_pockets[active_site_pocket_key]

    # step 4: identify possible fes pockets
    fes_pockets = identify_possible_fes_pockets(
        structure_db_path=structure_db_path,
        out=output_dir,
        chain_paths=chain_paths, 
        chain_dir=chain_dir_path,
        nn_clustering_radius=nn_clustering_radius,
        fident_threshold=fident_threshold,
        bits_threshold=bits_threshold
    )

    if fes_pockets == None:
        print("No iron-sulfur pockets found.")
        exit(0)

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
        if dist > nn_clustering_radius:
            valid_keys.append(pocket_key)
    for key in list(
        fes_pockets.keys()
    ):  # list() needed since you're modifying while iterating
        if key not in valid_keys:
            del fes_pockets[key]

    # step 6: check whether the pockets are really fes pockets using the fes pocket model
    fes_pocket_model, _, fingerprint_type, f_radius = load_model(fes_pocket_model_path)

    fes_pocket_probabilites = predict_cofactor_prob_for_pockets(
        pockets=fes_pockets,
        search_type=fes_pocket_search_type,
        structure=input_structure_path,
        model=fes_pocket_model,
        fingerprint_types=fingerprint_type,
        f_radius=f_radius
    )

    # step 7: filter fes pockets by pocket probability
    fes_pockets = _filter_fes_pockets_by_probability(fes_pockets, fes_pocket_probabilites)

    # step 8: predict each pocket
    fes_type_model, _, fingerprint_type, f_radius = load_model(fes_type_model_path)

    fes_probabilities = predict_cofactor_prob_for_pockets(
        pockets=fes_pockets,
        search_type=fes_type_search_type,
        structure=input_structure_path,
        model=fes_type_model,
        fingerprint_types=fingerprint_type,
        f_radius=f_radius
    )

    # step 9: 
    best_fes_hits_per_pocket = _choose_best_fes_per_pocket(
                    prediction_per_pocket=fes_probabilities
                )

    # step 10: fit sort and fit into hyd
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

    proximal_key, medial_key, distal_key = _fit_into_hyd_structure(
                    fes_cluster_distances=sorted_list,
                )

    write_outputs(
        out=output_dir,
        sorted_list=sorted_list,
        structure_path=input_structure_path,
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



##########################################################################################################################

def _fit_into_hyd_structure(fes_cluster_distances: List) -> tuple:
    """ """

    if len(fes_cluster_distances) == 3:
        proximal_key = fes_cluster_distances[0][0]
        medial_key = fes_cluster_distances[1][0]
        distal_key = fes_cluster_distances[2][0]
    else:
        # the proximal is around 10 +-
        proximal_key = None
        # the medial is around 20 +-
        medial_key = None
        # distal is around 30 +-
        distal_key = None
        for item in fes_cluster_distances:
            p = (12 - item[1])**2
            m = (20 - item[1])**2
            d = (30 - item[1])**2

            if p < m and p < d:
                if proximal_key != None:
                    print("Double proximal assignment")
                proximal_key = item[0]
            elif m < p and m < d:
                if medial_key != None:
                    print("Double medial assignment")
                medial_key = item[0]
            else:
                if distal_key != None:
                    print("Double distal assignment")
                distal_key = item[0]

    return proximal_key, medial_key, distal_key

def _choose_highest_conf_active_site_pocket(prediction_per_pocket: dict):
    """
    """
    best_hit_key = None
    best_hit_prob = 0.0

    for key in prediction_per_pocket:
        if best_hit_key == None:
            best_hit_key = key
            best_hit_prob = prediction_per_pocket[key]
            continue
        else:
            if best_hit_prob < prediction_per_pocket[key]:
                best_hit_key = key
                best_hit_prob = prediction_per_pocket[key]

    printl(f"Chose highest confidence active site with probability {best_hit_prob}")
    return best_hit_key

def _filter_fes_pockets_by_probability(fes_pockets, fes_probabilites):
    """"""
    res = {}
    for p in fes_probabilites.keys():
        if _get_best_hits_name(fes_probabilites[p]) == "fes_pocket":
            res[p] = fes_pockets[p]

    return res

def _get_best_hits_name(predictions: dict) -> str:
    """"""
    res = None
    res_logit = 0.0

    for k in predictions.keys():
        if predictions[k] > res_logit:
            res = k
            res_logit = predictions[k]

    return res

def _choose_best_fes_per_pocket(prediction_per_pocket: dict):
    """ """
    best_hits = {}

    for p in prediction_per_pocket.keys():
        best_hits[p] = _get_best_hits_name(prediction_per_pocket[p])

    return best_hits




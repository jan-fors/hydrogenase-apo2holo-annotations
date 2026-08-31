from typing import List
import numpy as np
from collections import defaultdict, Counter

from src.utils.geometric.calculate_geometric_centers import calculate_center
from src.fingerprint.create_fingerprint import create_fingerprint, create_fingerprints
from src.models.inference.predict import predict_pocket
from src.io.writers.printl import print_probabilities


def perform_mc_prediction(
    structure,
    model,
    fingerprint_types: List[str],
    f_radius: float,
    pocket_points: tuple[float, float, float],
):
    """ """
    M = np.array(pocket_points)
    center = calculate_center(M)
    F = [
        create_fingerprint(
            structure_path=structure,
            point=center,
            fingerprint_radius=f_radius,
            fingerprint_types=fingerprint_types,
        )
    ]
    pred = predict_pocket(model, F)
 
    return pred[0]


def perform_sum_prediction(
    structure,
    model,
    fingerprint_types: List[str],
    f_radius: float,
    pocket_points: tuple[float, float, float],
):
    """ """
    fingerprints = create_fingerprints(
        input_path=structure,
        coordinates=pocket_points,
        f_radius=f_radius,
        fingerprint_types=fingerprint_types,
    )
    pred = predict_pocket(model, fingerprints)

    return _aggregate_probs(pred)


def predict_cofactor_prob_for_pockets(
    structure,
    model,
    fingerprint_types: List[str],
    f_radius: float,
    search_type: str,
    pockets: dict,
):
    """ """
    pocket_results = {}

    # search each pocket
    for key in pockets.keys():
        pocket_results[key] = (
            perform_mc_prediction(
                structure, model, fingerprint_types, f_radius, pockets[key]
            )
            if search_type == "mc"
            else perform_sum_prediction(
                structure, model, fingerprint_types, f_radius, pockets[key]
            )
        )

        print("=" * 30, f"PROBS [{key}]", "=" * 30)
        print_probabilities(pocket_results[key])
        print("=" * 67)

    return pocket_results


def _aggregate_probs(prob_list):
    """ """
    acc = defaultdict(float)

    for p in prob_list:
        for k, v in p.items():
            acc[k] += v

    total = sum(acc.values())
    return {k: v / total for k, v in acc.items()}


def _counter_to_probs(counter):
    total = sum(counter.values())
    if total == 0:
        return {k: 0 for k in counter}
    return {k: v / total for k, v in counter.items()}
from pathlib import Path
from src.io.writers.printl import printl
from src.fingerprint.create_fingerprint import create_fingerprints
from src.models.inference.predict import predict_pocket

def predict_active_site_prob_by_pockets(
    pockets: dict,
    search_type: str,
    input_path: Path,
    model: Path = None,
    f_radius: float = None,
):
    """ """
    pocket_results = {}

    for key in pockets.keys():
        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = create_fingerprints(
                input_path=input_path, coordinates=pockets[key], f_radius=f_radius
            )

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = active_site_fingerprintDB.as_search(
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
            hits = active_site_fingerprintDB.as_search(
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
            hits = active_site_fingerprintDB.as_search(
                F=F, search_type=search_type_command, model=model
            )

            pocket_results[key] = _counter_to_probs(Counter(hits[0]))
        
        print("="*30, "PROBS", "="*30)
        print_probabilities(pocket_results[key])
        print("="*67)

    # return only the probability for an active site
    for key in pocket_results:
        pocket_results[key] = pocket_results[key]["active_site"]

    return pocket_results

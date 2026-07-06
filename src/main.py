import os
from pathlib import Path
from src.utils.protein.get_chains import get_chains
from src.utils.protein.extract_chain import extract_chain
from src.parser.parse_msearch_output import parse_msearch_output
from src.filter.filter_msearch_output import filter_msearch_output
from src.utils.protein.extract_cofactors import extract_cofactors
from src.utils.protein.apply_transformations_to_cofactors import apply_transformations_to_cofactors
import numpy as np
from src.parser.parse_t import parse_t
from src.parser.parse_u import parse_u
from src.utils.geometric.calculate_geometric_centers import calculate_geometric_centers, calculate_center
from src.io.plot import plot_protein
from src.filter.apply_blacklist import apply_blacklist
from src.filter.apply_whitelist import apply_whitelist, apply_active_site_whitelist, apply_fes_whitelist
from src.fingerprint.create_fingerprint import create_aminoacid_fingerprint, create_physiochemical_radial_angular, create_feature_fingerprint
from src.io.printl import printl, print_probabilities
from src.io.result_table import write_to_result_table
import json
from src.utils.geometric.nearest_neighbor_clustering import nn_radius_clustering
from tqdm import tqdm
from src.utils.protein.identfy_cysteines import identify_cysteines
from src.utils.geometric.point_distance import point_distance
from collections import Counter
from src.utils.constants import (
    STRUCTURE_DB,
    FINGERPRINT_DB,
    VERBOSE,
    SEARCH_TYPE,
    NN_CLUSTERING_RADIUS
)
from pprint import pprint
from src.database.StructureDB import StructureDB
from src.database.FingerprintDB import FingerprintDB
from src.utils.smiles import SMILES
from typing import List
from collections import defaultdict

create_fingerprint = create_physiochemical_radial_angular

######################## HELPER FUNCTIONS ########################

def _counter_to_probs(counter):
    total = sum(counter.values())
    if total == 0:
        return {k: 0 for k in counter}
    return {k: v / total for k, v in counter.items()}

def _aggregate_probs(prob_list):
    """
    """
    acc = defaultdict(float)

    for p in prob_list:
        for k, v in p.items():
            acc[k] += v

    total = sum(acc.values())
    return {k: v / total for k, v in acc.items()}

def _extract_protein_chains_from_file(input_path : Path, out : Path) -> List[Path]:
    """
    """
    # extract chains
    printl("Extract chains from protein file")
    chains = get_chains(input_path)
    chain_paths = []
    for chain in chains:
        chain_paths.append(extract_chain(input_path, out, chain))

    return chain_paths

def _structure_db_search(structureDB : StructureDB, chain_paths : List[Path], out : Path) -> List[tuple[float, float, float]]:
    """
    """
    # perform foldseek searches
    fd_res = []
    for cp in tqdm(chain_paths, disable=not VERBOSE):
        fd_res.append(structureDB.search(cp, out))

    for chain in chain_paths:
        # remove chain to save memory
        os.remove(chain)
        

    # parse the results
    cofactor_sites = []
    names = [] # ONLY relevant if plot is created

    printl("parse the results of the homology search in order to identify possible binding pockets")
    for fd_r in tqdm(fd_res, disable=not VERBOSE):
        # parse output
        r = parse_msearch_output(fd_r)

        # apply filters
        subset = filter_msearch_output(r)

        # extract coordinates for cofactors
        for i in range(subset.shape[0]):
            target = subset.iloc[i]["target"]
            u = subset.iloc[i]["u"]
            u_vec = parse_u(u)
   
            t = subset.iloc[i]["t"]
            t_vec = parse_t(t)

            # get structure path
            structure_path = structureDB.get_structure_path(target)
            if structure_path == None:
                continue

            # get cofactor coordinates
            cofactors = extract_cofactors(structure_path)

            cofactors = apply_blacklist(cofactors) #TODO change for active site search and fes search
            cofactors = apply_whitelist(cofactors)

            # apply transformations
            cofactors = apply_transformations_to_cofactors(cofactors, u_vec, t_vec)

            # calculate mass center
            geometric_centers = calculate_geometric_centers(cofactors)

            # get coordinate values
            coords = list(geometric_centers.values())
            names.extend(list(geometric_centers.keys()))
            cofactor_sites.extend(coords)

        # remove other results
        os.remove(fd_r)

    return cofactor_sites

def _active_site_structure_db_search(structureDB : StructureDB, chain_paths : List[Path], out : Path) -> List[tuple[float, float, float]]:
    """
    """
    
    # perform foldseek searches
    fd_res = []
    for cp in tqdm(chain_paths, disable=not VERBOSE):
        fd_res.append(structureDB.search(cp, out))

    # parse the results
    active_site_sites = []
    names = [] # ONLY relevant if plot is created

    printl("Parse the results of the active site homology search in order to identify possible binding pockets")
    for fd_r in tqdm(fd_res, disable=not VERBOSE):
        # parse output
        r = parse_msearch_output(fd_r)

        # apply filters
        subset = filter_msearch_output(r)

        # extract coordinates for cofactors
        for i in range(subset.shape[0]):
            target = subset.iloc[i]["target"]
            u = subset.iloc[i]["u"]
            u_vec = parse_u(u)
   
            t = subset.iloc[i]["t"]
            t_vec = parse_t(t)

            # get structure path
            structure_path = structureDB.get_structure_path(target)
            if structure_path == None:
                continue

            # get cofactor coordinates
            cofactors = extract_cofactors(structure_path)

            cofactors = apply_active_site_whitelist(cofactors)

            # apply transformations
            cofactors = apply_transformations_to_cofactors(cofactors, u_vec, t_vec)

            # calculate mass center
            geometric_centers = calculate_geometric_centers(cofactors)

            # get coordinate values
            coords = list(geometric_centers.values())
            names.extend(list(geometric_centers.keys()))
            active_site_sites.extend(coords)

        # remove other results
        os.remove(fd_r)

    return active_site_sites

def _fes_structure_db_search(structureDB : StructureDB, chain_paths : List[Path], out : Path) -> List[tuple[float, float, float]]:
    """
    """
    # perform foldseek searches
    fd_res = []
    for cp in tqdm(chain_paths, disable=not VERBOSE):
        fd_res.append(structureDB.search(cp, out))

    # parse the results
    fes_sites = []
    names = [] # ONLY relevant if plot is created

    printl("Parse the results of the active site homology search in order to identify possible binding pockets")
    for fd_r in tqdm(fd_res, disable=not VERBOSE):
        # parse output
        r = parse_msearch_output(fd_r)

        # apply filters
        subset = filter_msearch_output(r)

        # extract coordinates for cofactors
        for i in range(subset.shape[0]):
            target = subset.iloc[i]["target"]
            u = subset.iloc[i]["u"]
            u_vec = parse_u(u)
   
            t = subset.iloc[i]["t"]
            t_vec = parse_t(t)

            # get structure path
            structure_path = structureDB.get_structure_path(target)
            if structure_path == None:
                continue

            # get cofactor coordinates
            cofactors = extract_cofactors(structure_path)

            cofactors = apply_fes_whitelist(cofactors)

            # apply transformations
            cofactors = apply_transformations_to_cofactors(cofactors, u_vec, t_vec)

            # calculate mass center
            geometric_centers = calculate_geometric_centers(cofactors)

            # get coordinate values
            coords = list(geometric_centers.values())
            names.extend(list(geometric_centers.keys()))
            fes_sites.extend(coords)

        # remove other results
        os.remove(fd_r)

    for chain in chain_paths:
        # remove chain to save memory
        os.remove(chain)

    return fes_sites

def _arrange_cofactors_into_clusters(cofactor_sites : List[tuple[float, float, float]]) -> dict:
    """
    """
    # perform clustering
    labels = nn_radius_clustering(cofactor_sites, radius=NN_CLUSTERING_RADIUS)

    # sort into clusters
    cluster = {}
    for i in tqdm(range(len(labels)), disable=not VERBOSE):
        if labels[i] not in cluster.keys():
            cluster[labels[i]] = []
        cluster[labels[i]].append(cofactor_sites[i])

    return cluster

def _create_fingerprints(input_path : Path, coordinates : List[tuple[float, float, float]], f_radius : float = None):
    """"""
    fingerprints = []
    for point in tqdm(coordinates, disable=not VERBOSE):
        # create fingerprint
        if f_radius != None:
            F = create_fingerprint(input_path, point, f_radius)
        else:
            F = create_fingerprint(input_path, point)
        fingerprints.append(F)

    return fingerprints

def _get_best_hits_name(predictions : dict) -> str:
    """"""
    res = None
    res_logit = 0.0

    for k in predictions.keys():
        if predictions[k] > res_logit:
            res = k
            res_logit = predictions[k]

    return res

def _identify_active_site(best_hits : dict) -> int:
    """
    Edgecases:
    1. no active site -> Return None
    2. more than one acitve site -> Return -1
    """
    active_site_counter = 0
    active_site_index = None

    for h in best_hits:
        if best_hits[h] == "active_site":
            active_site_counter += 1
            active_site_index = int(h)

    if active_site_counter > 1:
        return -1
    else:
        return active_site_index
    
def _identify_fes_cluster(best_hits : dict) -> List[int]:
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

def _sort_into_hyd_format(active_site_mass_center : tuple, fes_cluster_mass_center : dict) -> tuple:
    """
    """
    proximal_key = None
    medial_key = None
    distal_key = None

    distances = []

    for k in fes_cluster_mass_center.keys():
        dist = point_distance(active_site_mass_center[0], active_site_mass_center[1], active_site_mass_center[2], fes_cluster_mass_center[k][0], fes_cluster_mass_center[k][1], fes_cluster_mass_center[k][2])
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

def identify_possible_cofactor_pockets(structure_db_path : Path, input_path : Path, out : Path):
    """
    The structure db is searched for good matches and the binding pockets of the template structures are printed onto the structure in question.
    Return:
        dict cluster:
            - list [cofactor sites x, y, z]
    """
    # init structure db
    printl("Initializing structure database...")
    structureDB = StructureDB()
    structureDB.load(structure_db_path=structure_db_path)

    # extract chains
    printl("Extracting protein chains...")
    chain_paths = _extract_protein_chains_from_file(input_path=input_path, out=out)

    # search against structuredb
    printl("Search against structure database for structural homologs...")
    cofactor_sites = _structure_db_search(structureDB=structureDB, chain_paths=chain_paths, out=out)

    # little check if something was found
    if len(cofactor_sites) == 0:
        print("No Cofactor Sites found")
        exit(0)
    else:
        printl(f"Found {len(cofactor_sites)} possible pockets")

    # sort cofactors into cluster through nearest neighbor clustering
    printl("Cluster cofactor sites into possible binding pockets...")
    clusters = _arrange_cofactors_into_clusters(cofactor_sites=cofactor_sites)
    printl(f"Sorted possible pockets into {len(clusters)} clusters.")

    return clusters

def predict_cofactors_by_pockets(fingerprint_db_path : Path, pockets : dict, search_type : str, input_path : Path, model : Path = None, f_radius : float = None):
    """
    """
    # Initialize fingerprintDB
    printl("Initializing fingerprint database...")
    fingerprintDB = FingerprintDB(fingerprint_db_path)

    
    # search each pocket against 
    printl("Searching each pocket point against the fingerprint database...")

    pocket_hits = {}

    for key in pockets.keys():

        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = fingerprintDB.search(F=fingerprints, search_type=search_type, model=model)

            # count
            pocket_hits[key] = Counter()
            pocket_hits[key].update(hits)

            pocket_hits[key] = _counter_to_probs(pocket_hits[key])

        elif search_type == "logreg_sum" or search_type == "mlp_sum" or search_type == "svm_sum" or search_type == "randforest_sum" or search_type == "sum":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fingerprintDB.search(F=fingerprints, search_type=search_type_command, model=model)

            ############# REMOVE
            AA_COLS_A = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
           'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
           'T', 'W', 'Y', 'V']

            print(AA_COLS_A)
            for f in fingerprints:
                print(f, f.sum())
            

            ############# REMOVE

            # count
            pocket_hits[key] = _aggregate_probs(hits)

        elif search_type == "logreg_mc" or search_type == "mlp_mc" or search_type == "svm_mc" or search_type == "randforest_mc" or search_type == "mc":
            # create center
            M = np.array(pockets[key])
            center = calculate_center(M)
            printl(f"Creating fingerprint for pocket {key}...")
            if f_radius != None:
                F = [create_fingerprint(input_path, center, f_radius)]
            else:
                F = [create_fingerprint(input_path, center)]

            ############# REMOVE
            AA_COLS_A = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
           'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
           'T', 'W', 'Y', 'V']

            print(AA_COLS_A)
            print(F, np.array(F).sum())
            

            ############# REMOVE
            
            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fingerprintDB.search(F=F, search_type=search_type_command, model=model)

            pocket_hits[key] = _counter_to_probs(Counter(hits[0]))



        print_probabilities(pocket_hits[key])

    return pocket_hits

def choose_best_hit_per_pocket(prediction_per_pocket : dict) -> dict:
    """
    """
    best_hits = {}

    for p in prediction_per_pocket.keys():
        best_hits[p] = _get_best_hits_name(prediction_per_pocket[p])

    return best_hits

def fit_into_hyd_structure(active_site_mass_center : tuple, distances : List, best_hits : dict, pockets : dict) -> tuple:
    """
    """
    

    fes_cluster_keys = _identify_fes_cluster(best_hits)
    fes_cluster_mass_centers = {}
    for k in fes_cluster_keys:
        fes_cluster_mass_centers[k] = calculate_center(np.array(pockets[k]))    
    
    proximal_key, medial_key, distal_key = _sort_into_hyd_format(active_site_mass_center, fes_cluster_mass_centers)

    return proximal_key, medial_key, distal_key

def write_outputs(out : Path, sorted_list : List, structure_path : Path, as_pockets : dict, as_pocket_key : int, as_probabilities : dict, fes_pockets : dict,                   prediction_per_pocket : dict, best_hits_per_pocket : dict, proximal_key : int, medial_key : int, distal_key : int, boltz : bool, plot : bool):
    """"""
    # raw_results
    output_dict = {}

    # active_site
    output_dict["active_site"] = {
        "prediction": as_probabilities[as_pocket_key],
        "cysteines": identify_cysteines(structure_path=structure_path, coords=calculate_center(np.array(as_pockets[as_pocket_key])), radius=5.0)
    }

    for k in fes_pockets.keys():
        #print(k)
        if best_hits_per_pocket[k] != "protein":
            cysteines = identify_cysteines(structure_path=structure_path, coords=calculate_center(np.array(fes_pockets[k])), radius=5.0)
            output_dict[k] = {
                "cluster": int(k),
                "predicted_class": best_hits_per_pocket[k],
                "smiles": SMILES[best_hits_per_pocket[k]],
                "prediction": prediction_per_pocket[k],
                "cysteines": cysteines
            }


    output_list = []
    output_list.append(output_dict["active_site"])

    for i in sorted_list:
        if int(i[0]) in output_dict.keys():
            output_list.append(output_dict[i[0]])

    ## write json output file
    json_path = os.path.join(out, "raw_results.json")
    json_text = json.dumps(output_list, ensure_ascii=False, indent=4)
    printl(f"Writing raw output to {json_path}")
    with open(json_path, "w", encoding="utf-8") as f:
        f.write(json_text)


    # results
    result_dict = {}

    result_dict["active_site"] = output_dict["active_site"],
    if proximal_key != None:
        result_dict["proximal"] = output_dict[proximal_key],
    if medial_key != None:
        result_dict["medial"] = output_dict[medial_key],
    if distal_key != None:
        result_dict["distal"] = output_dict[distal_key]
    
    if len(result_dict.keys()) > 0:
        json_path = os.path.join(out, "results.json")
        json_text = json.dumps(result_dict, ensure_ascii=False, indent=4)
        printl(f"Writing structured output to {json_path}")
        with open(json_path, "w", encoding="utf-8") as f:
            f.write(json_text)

    # write boltz yaml output file
    if boltz:
        print("YAML CREATION TODO")

    # create output graphic?
    if plot:
        plot_protein(structure_path = structure_path, pockets = fes_pockets, pred_per_pocket = prediction_per_pocket, output = out)

def _identify_possible_active_site_pockets(structure_db_path : Path, input_path : Path, out : Path) -> tuple[list, dict]:
    """
    Identifies the active_site pockets inside the input_structure
    Returns the chain paths and a dictionary, where the keys are the clusters and the items are the pocket points
    """

    # init active-site structure db
    printl("Initializing active-site structure database ...")
    active_site_structureDB = StructureDB()
    active_site_structure_db_path = structure_db_path
    active_site_structureDB.load(structure_db_path=active_site_structure_db_path)

    # extract chains
    printl("Extracting protein chains ...")
    chain_paths = _extract_protein_chains_from_file(input_path=input_path, out=out)

    # search against active-site structure db
    printl("Search against active-site structure database for structural homologs ...")
    active_site_sites = _active_site_structure_db_search(structureDB=active_site_structureDB, chain_paths=chain_paths, out=out)

    # intermediate check
    if len(active_site_sites) == 0:
        print("No possible active sites found")
        exit(0)
    else:
        printl(f"Found {len(active_site_sites)} possible active-site pockets")

    # sort into clusters
    printl("Cluster cofactor sites into possible binding pockets ...")
    clusters = _arrange_cofactors_into_clusters(cofactor_sites=active_site_sites)
    printl(f"Sorted possible pockets into {len(clusters)} clusters.")

    return chain_paths, clusters

def _predict_active_site_prob_by_pockets(fingerprint_db_path : Path, pockets : dict, search_type : str, input_path : Path, model : Path = None, f_radius : float = None):
    """
    """
    # Initialize fingerprintDB
    printl("Initializing active-site fingerprint database...")
    active_site_fingerprint_db_path = fingerprint_db_path / Path("activesiteDB")
    active_site_fingerprintDB = FingerprintDB(active_site_fingerprint_db_path)

    # search each pocket against 
    printl("Searching each pocket against the active-site fingerprint database ...")

    pocket_results = {}

    for key in pockets.keys():
        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = active_site_fingerprintDB.as_search(F=fingerprints, search_type=search_type, model=model)

            # count
            pocket_results[key] = Counter()
            pocket_results[key].update(hits)

            pocket_results[key] = _counter_to_probs(pocket_results[key])

        elif search_type == "logreg_sum" or search_type == "mlp_sum" or search_type == "svm_sum" or search_type == "randforest_sum" or search_type == "sum":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = active_site_fingerprintDB.as_search(F=fingerprints, search_type=search_type_command, model=model)

            # count
            pocket_results[key] = _aggregate_probs(hits)

        elif search_type == "logreg_mc" or search_type == "mlp_mc" or search_type == "svm_mc" or search_type == "randforest_mc" or search_type == "mc":
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
            hits = active_site_fingerprintDB.as_search(F=F, search_type=search_type_command, model=model)

            pocket_results[key] = _counter_to_probs(Counter(hits[0]))

        print_probabilities(pocket_results[key])

    # return only the probability for an active site
    for key in pocket_results:
        pocket_results[key] = pocket_results[key]["active_site"]

    return pocket_results

def _choose_highest_conf_active_site_pocket(prediction_per_pocket : dict):
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

def _identify_possible_fes_pockets(structure_db_path : Path, chain_paths : List, input_path : Path, out : Path):
    """
    """
    # initialize database
    printl("Initialize fes structure database ...")
    fes_structureDB = StructureDB()
    fes_structure_db_path = structure_db_path
    fes_structureDB.load(structure_db_path=fes_structure_db_path)

    # search against structure db
    printl("Search against fes structure database for structural homologs ...")
    fes_sites = _fes_structure_db_search(structureDB=fes_structureDB, chain_paths=chain_paths, out=out)

    # intermediate check
    if len(fes_sites) == 0:
        print("No possible active sites found")
        exit(0)
    else:
        printl(f"Found {len(fes_sites)} possible fes pockets")

    # sort into clusters
    printl("Cluster cofactor sites into possible binding pockets ...")
    clusters = _arrange_cofactors_into_clusters(cofactor_sites=fes_sites)
    printl(f"Sorted possible pockets into {len(clusters)} clusters.")

    return clusters

def _predict_fes_prob_by_pockets(fingerprint_db_path : Path, pockets : dict, search_type : str, input_path : Path, model : Path = None, f_radius : float = None):
    # Initialize fingerprintDB
    printl("Initializing fes fingerprint database...")
    fes_fingerprint_db_path = fingerprint_db_path / Path("fingerprintDB") #TODO is this needed?
    fes_fingerprintDB = FingerprintDB(fes_fingerprint_db_path)

    # search each pocket against 
    printl("Searching each pocket against the active-site fingerprint database ...")

    pocket_results = {}

    for key in pockets.keys():
        if search_type == "absolut":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")
            hits = fes_fingerprintDB.fes_search(F=fingerprints, search_type=search_type, model=model)

            # count
            pocket_results[key] = Counter()
            pocket_results[key].update(hits)

            pocket_results[key] = _counter_to_probs(pocket_results[key])

        elif search_type == "logreg_sum" or search_type == "mlp_sum" or search_type == "svm_sum" or search_type == "randforest_sum" or search_type == "sum":
            printl(f"Creating fingerprints for pocket {key}...")
            fingerprints = _create_fingerprints(input_path=input_path, coordinates=pockets[key], f_radius=f_radius)

            # search
            printl(f"Searching pocket {key} against fingerprint database...")

            if "_" in search_type:
                search_type_command = search_type.split("_")[0]
            else:
                search_type_command = search_type
            hits = fes_fingerprintDB.fes_search(F=fingerprints, search_type=search_type_command, model=model)

            # count
            pocket_results[key] = _aggregate_probs(hits)

        elif search_type == "logreg_mc" or search_type == "mlp_mc" or search_type == "svm_mc" or search_type == "randforest_mc" or search_type == "mc":
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
            hits = fes_fingerprintDB.fes_search(F=F, search_type=search_type_command, model=model)

            pocket_results[key] = _counter_to_probs(Counter(hits[0]))

        print_probabilities(pocket_results[key])

    return pocket_results

def _choose_best_fes_per_pocket(prediction_per_pocket : dict):
    """
    """
    best_hits = {}

    for p in prediction_per_pocket.keys():
        best_hits[p] = _get_best_hits_name(prediction_per_pocket[p])

    return best_hits

def main(input_path : Path, out : Path, tmp : Path, boltz : bool, plot: bool, structure_db_path : Path = STRUCTURE_DB, fingerprint_db_path : Path = FINGERPRINT_DB, search_type=SEARCH_TYPE, active_site_model : Path = None, fes_model : Path = None, f_radius : float = None):
    """

    """
    # step 1: identify possible active site pockets
    chain_paths, active_site_pockets = _identify_possible_active_site_pockets(structure_db_path=structure_db_path, input_path=input_path, out=out)

    # step 2: predict active site probability per pocket
    active_site_probabilites = _predict_active_site_prob_by_pockets(fingerprint_db_path=fingerprint_db_path, pockets=active_site_pockets, search_type=search_type, input_path=input_path, model=active_site_model, f_radius=f_radius)

    # step 3: choose highest confidence active site
    active_site_pocket_key = _choose_highest_conf_active_site_pocket(active_site_probabilites)
    active_site_cluster_points = active_site_pockets[active_site_pocket_key]

    # step 4: identify possible fes pockets
    fes_pockets = _identify_possible_fes_pockets(structure_db_path=structure_db_path, chain_paths=chain_paths, input_path=input_path, out=out)

    # step 5: check if clusters overlap with active_site cluster
    active_site_cluster_center = calculate_center(np.array(active_site_cluster_points))

    valid_keys = []
    for pocket_key in fes_pockets.keys():
        pocket_center = calculate_center(np.array(fes_pockets[pocket_key]))
        
        dist = point_distance(active_site_cluster_center[0], active_site_cluster_center[1], active_site_cluster_center[2], pocket_center[0], pocket_center[1], pocket_center[2])
        if dist > NN_CLUSTERING_RADIUS:
            valid_keys.append(pocket_key)
    for key in list(fes_pockets.keys()):  # list() needed since you're modifying while iterating
        if key not in valid_keys:
            del fes_pockets[key]

    # step 6: predict each pocket
    fes_probabilities = _predict_fes_prob_by_pockets(fingerprint_db_path=fingerprint_db_path, pockets=fes_pockets, search_type=search_type, input_path=input_path, model=fes_model, f_radius=f_radius)
    
    # step 7: choose highest confidence hits
    best_fes_hits_per_pocket = _choose_best_fes_per_pocket(prediction_per_pocket=fes_probabilities)

    # step 8: fit sort and fit into hyd
    cluster_dist = []
    for pocket_key in fes_pockets.keys():
        pocket_center = calculate_center(np.array(fes_pockets[pocket_key]))
        dist = point_distance(active_site_cluster_center[0], active_site_cluster_center[1], active_site_cluster_center[2], pocket_center[0], pocket_center[1], pocket_center[2])
        cluster_dist.append((pocket_key, dist))

    sorted_list = sorted(cluster_dist, key=lambda t: t[1])

    proximal_key, medial_key, distal_key = fit_into_hyd_structure(active_site_mass_center=active_site_cluster_center ,distances=sorted_list, best_hits=best_fes_hits_per_pocket, pockets=fes_pockets)
    
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
        plot=plot
        )

    



import os
from pathlib import Path
#from src.database.search_against_structuredb import search_against_structuredb
from src.utils.get_chains import get_chains
from src.utils.extract_chain import extract_chain
from src.parser.parse_msearch_output import parse_msearch_output
from src.filter.filter_msearch_output import filter_msearch_output
#from src.database.get_structure_path import get_structure_path
from src.utils.extract_cofactors import extract_cofactors
from src.utils.apply_transformations_to_cofactors import apply_transformations_to_cofactors
import numpy as np
from src.parser.parse_t import parse_t
from src.parser.parse_u import parse_u
from src.utils.calculate_geometric_centers import calculate_geometric_centers, calculate_center
from src.io.plot import plot_with_protein_from_pdb
from src.filter.apply_blacklist import apply_blacklist
from src.filter.apply_whitelist import apply_whitelist
from src.fingerprint.create_fingerprint import create_fingerprint
from src.io.printl import printl
from src.io.result_table import write_to_result_table
import json
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components
from tqdm import tqdm
from src.utils.identfy_cysteines import identify_cysteines
from src.utils.point_distance import point_distance
from collections import Counter
from src.utils.constants import (
    STRUCTURE_DB,
    FINGERPRINT_DB,
    VERBOSE,
    SEARCH_TYPE
)
from pprint import pprint
from src.database.StructureDB import StructureDB
from src.database.FingerprintDB import FingerprintDB
from src.utils.smiles import SMILES
from collections import defaultdict

def nn_radius_clustering(points, radius): #TODO move to own file
    """
    points: numpy array (N,3)
    radius: distance threshold

    Returns:
        labels: numpy array (N,)
    """

    # Build radius neighbor graph
    nbrs = NearestNeighbors(radius=radius)
    nbrs.fit(points)

    adjacency_matrix = nbrs.radius_neighbors_graph(points)

    # Find connected components
    n_components, labels = connected_components(adjacency_matrix)

    return labels

def aggregate_probs(prob_list):
    acc = defaultdict(float)

    for p in prob_list:
        for k, v in p.items():
            acc[k] += v

    total = sum(acc.values())
    return {k: v / total for k, v in acc.items()}

def main(input_path : str, output : str, output_dir : str, tmp : str, boltz : bool, structure_db_path : str = STRUCTURE_DB, fingerprint_db_path : str = FINGERPRINT_DB, search_type=SEARCH_TYPE, result_table_path : str = None):
    """
    1. Amino-acid sequence of structure is BLASTed against the sequence file of the database, which contains known hydrogenase structures and proteins that contain FeS-Cofactors. Return hits of database sorted by E values.
    2. Structurally align the hits with the input protein on the Ca-atoms of the residues that match in the BLAST alignment.
    3. Identify the positions of the Active Site as well as the FeS Cofactors
    4. Check the surrounding of the mass center and create aminoacid fingerprints
    5. Group the cofactors into the groups: activesite, proximal, medial and distal cluster
    6. select the most reasonable combination of cofactors based on the fingerprints.
    7. create output file and return
    """
    # create outputfolder with name output 
    # TODO move to cli.py
    out = Path(os.path.join(output_dir, output))
    os.makedirs(out, exist_ok=True)
    
    # extract chains
    printl("Extract chains from protein file")
    chains = get_chains(input_path)
    chain_paths = []
    for chain in chains:
        chain_paths.append(extract_chain(input_path, out, chain))
    
    # run foldseek against structure db
    printl("Initializing structure db")
    structureDB = StructureDB()
    structureDB.load(structure_db_path=structure_db_path)
    
    printl("Search chains against structure db for structural homologs")
    fd_res = []
    for cp in tqdm(chain_paths, disable=not VERBOSE):
        fd_res.append(structureDB.search(cp, out))
        #fd_res.append(search_against_structuredb(cp, out, tmp, structure_db_path))

    cofactor_sites = []
    names = []

    # parse results for each chain and extract cofactors
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

            cofactors = apply_blacklist(cofactors)
            cofactors = apply_whitelist(cofactors)

            # apply transformations
            cofactors = apply_transformations_to_cofactors(cofactors, u_vec, t_vec)

            # calculate mass center
            geometric_centers = calculate_geometric_centers(cofactors)

            # get coordinate valus
            coords = list(geometric_centers.values())
            names.extend(list(geometric_centers.keys()))
            cofactor_sites.extend(coords)

    if len(cofactor_sites) == 0:
        print("No Cofactor Sites found")
        exit(0)

    # sort cofactors into clusters
    printl("cluster the possible binding pockets")
    labels = nn_radius_clustering(cofactor_sites, radius=4.0)

    # rearrange data
    printl("Rearrange data into pockets")
    cluster = {}
    for i in tqdm(range(len(labels)), disable=not VERBOSE):
        if labels[i] not in cluster.keys():
            cluster[labels[i]] = []
        cluster[labels[i]].append(cofactor_sites[i])
        


    # fig = plot_with_protein_from_pdb(
    # pdb_path=input_path,
    # cofactor_coords=cofactor_sites,
    # labels=labels,
    # names=names,
    # out_html="protein_plot.html"
    # )

    # for each label -> search each point against fingerprintDB
    printl("Initialize Fingerprint DB")
    fingerprintDB = FingerprintDB(fingerprint_db_path)

    printl("Search each pocket point against the fingerprint db")
    results = {}
    for cl in cluster.keys():
        if SEARCH_TYPE == "absolut":
            results[cl] = Counter()

            # for each cofactor create fingerprints
            fingerprints = []
            for point in tqdm(cluster[cl], disable=not VERBOSE):
                # create fingerprint
                F = create_fingerprint(input_path, point)
                fingerprints.append(F)
                
            hits = fingerprintDB.search(fingerprints, SEARCH_TYPE)
            # count occurrences
            results[cl].update(hits)

        elif SEARCH_TYPE == "logreg_mc":
            # calculate mass center of cluster
            M = np.array(cluster[cl])
            center = calculate_center(M)
            F = [create_fingerprint(input_path, center)]
            hits = fingerprintDB.search(F, "logreg")

            results[cl] = Counter(hits[0])

        elif SEARCH_TYPE == "logreg_sum":
            results[cl] = Counter()

            # for each cofactor create fingerprints
            fingerprints = []
            for point in tqdm(cluster[cl], disable=not VERBOSE):
                # create fingerprint
                F = create_fingerprint(input_path, point)
                fingerprints.append(F)

            # check fingerprintdb
            hits = fingerprintDB.search(fingerprints, "logreg")

            # count occurrences
            results[cl] = aggregate_probs(hits)

        elif SEARCH_TYPE == "svm_sum":
            results[cl] = Counter()

            # for each cofactor create fingerprints
            fingerprints = []
            for point in tqdm(cluster[cl], disable=not VERBOSE):
                # create fingerprint
                F = create_fingerprint(input_path, point)
                fingerprints.append(F)

            # check fingerprintdb
            hits = fingerprintDB.search(fingerprints, "svm")

            # count occurrences
            results[cl] = aggregate_probs(hits)

        elif SEARCH_TYPE == "svm_mc":
            # calculate mass center of cluster
            M = np.array(cluster[cl])
            center = calculate_center(M)
            F = [create_fingerprint(input_path, center)]
            hits = fingerprintDB.search(F, "svm")

            results[cl] = Counter(hits[0])

        elif SEARCH_TYPE == "mlp_sum":
            results[cl] = Counter()

            # for each cofactor create fingerprints
            fingerprints = []
            for point in tqdm(cluster[cl], disable=not VERBOSE):
                # create fingerprint
                F = create_fingerprint(input_path, point)
                fingerprints.append(F)

            # check fingerprintdb
            hits = fingerprintDB.search(fingerprints, "mlp")

            # count occurrences
            results[cl] = aggregate_probs(hits)

        elif SEARCH_TYPE == "mlp_mc":
            # calculate mass center of cluster
            M = np.array(cluster[cl])
            center = calculate_center(M)
            F = [create_fingerprint(input_path, center)]
            hits = fingerprintDB.search(F, "mlp")

            results[cl] = Counter(hits[0])

    print(results)
    exit(0)

    """
    The selection mechanic might change using a different database/search engine.
    For test reasons:
    - for every cluster that has hits
    - choose type that has most hits
        - FeS or Active Site
            - if active site cluster -> active site smiles
            - if FeS wins -> type smiles
    """

    # select for each cluster the best result
    printl("Select the best matching cofactors")

    active_site = {
        "x": None,
        "y": None,
        "z": None,
        "cysteines": None
    }

    # identify active site
    key_active_site = None
    amount_active_site = 0
    proba_active_site = 0.0
    for cl in results.keys():
        if SEARCH_TYPE == "absolut":
            tmp = dict(results[cl])
            if not tmp:
                continue

            lst = sorted(tmp.items(), key=lambda x: x[1], reverse=True)

            first = lst[0][0]
            try:
                second = lst[1][0]
            except:
                second = ""

            if "active_site" in first or "active_site" in second:
                if key_active_site == None:
                    key_active_site = cl
                    amount_active_site = tmp["active_site"]
                else:
                    if amount_active_site < tmp["active_site"]:
                        key_active_site = cl
                        amount_active_site = tmp["active_site"]

        elif SEARCH_TYPE == "logreg_mc":
            # get the one with the highest probability
            # check if results[cl] is empty
            if not results[cl]:
                continue
            if results[cl][0][0] == "active_site":
                if key_active_site == None:
                    key_active_site = cl
                    proba_active_site = results[cl][0][1]
                else:
                    if proba_active_site < results[cl][0][1]:
                        key_active_site = cl
                        proba_active_site = results[cl][0][1]

        elif SEARCH_TYPE == "logreg_sum":
            tmp = dict(results[cl])
            if not tmp:
                continue

            lst = sorted(tmp.items(), key=lambda x: x[1], reverse=True)

            first = lst[0][0]
            try:
                second = lst[1][0]
            except:
                second = ""

            if "active_site" in first or "active_site" in second:
                if key_active_site == None:
                    key_active_site = cl
                    amount_active_site = tmp["active_site"]
                else:
                    if amount_active_site < tmp["active_site"]:
                        key_active_site = cl
                        amount_active_site = tmp["active_site"]

    # assign active site
    active_site["formula"] = "active_site"
    # calculate cluster mass centers for clusters with  
    M = np.array(cluster[key_active_site])

    center = calculate_center(M)
    active_site["x"] = center[0]
    active_site["y"] = center[1]
    active_site["z"] = center[2]
    # identfy cysteins which are important for binding
    cysteines = identify_cysteines(input_path, center)
    active_site["cysteines"] = cysteines

    # go through other cluster
    fes_cluster = []
    for cl in results.keys():
        if cl == key_active_site:
            continue

        if SEARCH_TYPE == "absolut":
            tmp = dict(results[cl])
            if not tmp:
                print(f"{cl} is an empty cluster. Removing it..")
                continue

            lst = sorted(tmp.items(), key=lambda x: x[1], reverse=True) #TODO muss true sein?
            formula = lst[0][0]

            if formula == "active_site":
                formula = lst[1][0]
            
        elif SEARCH_TYPE == "logreg_mc":
            formula = results[cl][0][0]

        elif SEARCH_TYPE == "logreg_sum":
            tmp = dict(results[cl])
            if not tmp:
                print(f"{cl} is an empty cluster. Removing it..")
                continue

            lst = sorted(tmp.items(), key=lambda x: x[1], reverse=True) #TODO muss true sein?
            formula = lst[0][0]

            if formula == "active_site":
                formula = lst[1][0]

        # calculate cluster mass centers for clusters with  
        M = np.array(cluster[cl])

        center = calculate_center(M)

        # identfy cysteins which are important for binding
        cysteines = identify_cysteines(input_path, center)

        printl(f"Cluster: {cl}, formula: {formula}, cysteines: {cysteines}")

    
        if formula != "protein":
            fes_cluster.append(
                {
                    "formula": formula,
                    "x": center[0],
                    "y": center[1],
                    "z": center[2],
                    "cysteines": cysteines
                }
            )

    # try to map into known system
    is_conform = False

    # calculate distanz of each fes cluster to active site
    for i in range(len(fes_cluster)):
        fes_cluster[i]["dist"] = point_distance(active_site["x"], active_site["y"], active_site["z"], fes_cluster[i]["x"], fes_cluster[i]["y"], fes_cluster[i]["z"])

    # sort by length
    sorted_data = sorted(fes_cluster, key=lambda d: d["dist"])

    # validate/interpret output
    if len(fes_cluster) > 3:
        proximal = None
        medial = None
        distal = None
        printl("Identified more than three fes cluster. Cannot choose proximal, medial and distal.")
    elif len(fes_cluster) < 3:
        proximal = None
        medial = None
        distal = None
        printl("Identified less than 3 fes clusters. Cannot choose proximal, medial and distal.")
    else:
        is_conform = True
        
        proximal = sorted_data[0]
        medial = sorted_data[1]
        distal = sorted_data[2]

    # define result

    if is_conform:
        result = {
            "active_site": {
                "smiles": ":)",
                "formula": active_site["formula"],
                "coords":  {
                    "x": active_site["x"],
                    "y": active_site["y"],
                    "z": active_site["z"]
                    },
                "cystein-connections": active_site["cysteines"]
            },
            "proximal": {
                "smiles": ":)",
                "formula": proximal["formula"],
                "coords":  {
                    "x": proximal["x"],
                    "y": proximal["y"],
                    "z": proximal["z"]
                    },
                "cystein-connections": proximal["cysteines"]
            },
            "medial": {
                "smiles": ":)",

                "formula": medial["formula"],
                "coords":  {
                    "x": medial["x"],
                    "y": medial["y"],
                    "z": medial["z"]
                    },
                "cystein-connections": medial["cysteines"]
            },
            "distal": {
                "smiles": ":)",


                "formula": distal["formula"],
                "coords":  {
                    "x": distal["x"],
                    "y": distal["y"],
                    "z": distal["z"]
                    },
                "cystein-connections": distal["cysteines"]
            }
        }
    else:
        result = {
            "active_site": {
                "smiles": ":)",
                "formula": active_site["formula"],
                "coords":  {
                    "x": active_site["x"],
                    "y": active_site["y"],
                    "z": active_site["z"]
                    },
                "cystein-connections": active_site["cysteines"]
            }}
        for i in range(len(fes_cluster)):
            result[i] = {
                "smiles": ":)",
                "formula": fes_cluster[i]["formula"],
                "coords":  {
                    "x": fes_cluster[i]["x"],
                    "y": fes_cluster[i]["y"],
                    "z": fes_cluster[i]["z"]
                    },
                "cystein-connections": fes_cluster[i]["cysteines"]
            }

    pprint(result)

    # write json output file
    json_path = os.path.join(out, "result.json")
    json_text = json.dumps(result, ensure_ascii=False, indent=4)
    with open(json_path, "w", encoding="utf-8") as f:
        f.write(json_text)

    #   for testing write to table?
    if result_table_path:
        write_to_result_table(result, result_table_path)

    # create yaml TODO

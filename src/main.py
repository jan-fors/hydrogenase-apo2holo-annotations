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
from src.database.search_against_fingerprint_db import search_against_fingerprint_db
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
from src.database.StructureDB import StructureDB
from src.database.FingerprintDB import FingerprintDB
from src.utils.smiles import SMILES

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

def main(input_path : str, output : str, output_dir : str, tmp : str, boltz : bool, structure_db_path : str = STRUCTURE_DB, fingerprint_db_path : str = FINGERPRINT_DB, result_table_path : str = None):
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
    structureDB.load(structure_db_path=STRUCTURE_DB)
    
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
        


    fig = plot_with_protein_from_pdb(
    pdb_path=input_path,
    cofactor_coords=cofactor_sites,
    labels=labels,
    names=names,
    out_html="protein_plot.html"
    )

    # for each label -> search each point against fingerprintDB
    printl("Initialize Fingerprint DB")
    fingerprintDB = FingerprintDB(FINGERPRINT_DB)

    printl("Search each pocket point against the fingerprint db")
    results = {}
    for cl in cluster.keys():
        results[cl] = Counter()

        # for each cofactor
        for point in tqdm(cluster[cl], disable=not VERBOSE):
            # create fingerprint
            F = create_fingerprint(input_path, point)

            # check fingerprintdb
            #hits = search_against_fingerprint_db(F, fingerprint_db_path)
            if SEARCH_TYPE == "absolut":
                hits = fingerprintDB.search(F, SEARCH_TYPE)
            elif SEARCH_TYPE == "logreg":
                hits = fingerprintDB.search(F, SEARCH_TYPE)

                if hits:
                    hits = fingerprintDB.get_model().classes_[hits[0].argmax(axis=1)]
            # count occurrences
            results[cl].update(hits)

    print(results)
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

    fes_cluster = []
    for cl in results.keys():
        tmp = dict(results[cl])
        if not tmp:
            print(f"{cl} is an empty cluster. Removing it..")
            continue

        lst = sorted(tmp.items(), key=lambda x: x[1], reverse=True) #TODO muss true sein?
        formula = lst[0][0]

        try:
            formula_sec = lst[1][0] #TODO just quickfix
        except:
            formula_sec = ""

        is_active_site = False

        if "active_site" in formula or "active_site" in formula_sec:# in ["NFV", "NI", "NFU", "3NI", "FCO", "CMO"]: # active site
            is_active_site = True
            formula = "active_site"
     
        # calculate cluster mass centers for clusters with  
        M = np.array(cluster[cl])

        center = calculate_center(M)

        # identfy cysteins which are important for binding
        cysteines = identify_cysteines(input_path, center)

        printl(f"Cluster: {cl}, is active stie: {is_active_site}, formula: {formula}, cysteines: {cysteines}")

        if is_active_site:
            active_site["formula"] = formula
            active_site["x"] = center[0]
            active_site["y"] = center[1]
            active_site["z"] = center[2]
            active_site["cysteines"] = cysteines
        else:
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

    # validate/interpret output
    if len(fes_cluster) > 3:
        proximal = None
        medial = None
        distal = None
        printl("Identified more than three fes cluster. Cannot choose proximal, medial and distal.")
        exit(0)
    elif len(fes_cluster) < 3:
        proximal = None
        medial = None
        distal = None
        printl("Identified less than 3 fes clusters. Cannot choose proximal, medial and distal.")
        exit(0)
    else:
        # calculate distanz of each fes cluster to active site
        for i in range(len(fes_cluster)):
            fes_cluster[i]["dist"] = point_distance(active_site["x"], active_site["y"], active_site["z"], fes_cluster[i]["x"], fes_cluster[i]["y"], fes_cluster[i]["z"])

        # sort by length
        sorted_data = sorted(fes_cluster, key=lambda d: d["dist"])

        proximal = sorted_data[0]
        medial = sorted_data[1]
        distal = sorted_data[2]

    # define result


    result = {
        "active_site": {
            "smiles": ":)",
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

    # write json output file
    json_path = os.path.join(out, "result.json")
    json_text = json.dumps(result, ensure_ascii=False, indent=4)
    with open(json_path, "w", encoding="utf-8") as f:
        f.write(json_text)

    #   for testing write to table?
    if result_table_path:
        write_to_result_table(result, result_table_path)

    # create yaml TODO

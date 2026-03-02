import os
from pathlib import Path
from src.db.search_against_structuredb import search_against_structuredb
from src.utils.get_chains import get_chains
from src.utils.extract_chain import extract_chain
from src.parser.parse_msearch_output import parse_msearch_output
from src.filter.filter_msearch_output import filter_msearch_output
from src.db.get_structure_path import get_structure_path
from src.utils.extract_cofactors import extract_cofactors
from src.utils.apply_transformations_to_cofactors import apply_transformations_to_cofactors
import numpy as np
from src.parser.parse_t import parse_t
from src.parser.parse_u import parse_u
from src.utils.calculate_geometric_centers import calculate_geometric_centers
from src.io.plot import plot_with_protein_from_pdb
from src.filter.apply_blacklist import apply_blacklist

from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components


from sklearn.mixture import GaussianMixture

def nn_radius_clustering(points, radius):
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

def main(input_path : str, output : str, output_dir : str, tmp : str, boltz : bool):
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
    out = Path(os.path.join(output_dir, output))
    os.makedirs(out, exist_ok=True)
    
    # extract chains
    chains = get_chains(input_path)
    chain_paths = []
    for chain in chains:
        chain_paths.append(extract_chain(input_path, out, chain))
    
    # run foldseek against structure db
    fd_res = []
    for cp in chain_paths:
        fd_res.append(search_against_structuredb(cp, out, tmp))

    cofactor_sites = []
    names = []

    # parse results for each chain
    for fd_r in fd_res:
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
            structure_path = get_structure_path(target)
            
            # get cofactor coordinates
            cofactors = extract_cofactors(structure_path)

            cofactors = apply_blacklist(cofactors)

            # apply transformations
            cofactors = apply_transformations_to_cofactors(cofactors, u_vec, t_vec)

            # calculate mass center
            geometric_centers = calculate_geometric_centers(cofactors)

            # get coordinate valus
            coords = list(geometric_centers.values())
            names.extend(list(geometric_centers.keys()))
            cofactor_sites.extend(coords)

    #group coordinates
    # kmeans = KMeans(n_clusters=4, n_init="auto", random_state=0)
    # labels = kmeans.fit_predict(cofactor_sites)

    # gmm = GaussianMixture(n_components=4, covariance_type="full")
    # labels = gmm.fit_predict(cofactor_sites)

    # db = DBSCAN(eps=3.0, min_samples=5)
    # labels = db.fit_predict(cofactor_sites)

    labels = nn_radius_clustering(cofactor_sites, radius=5.0)

    for i in range(len(labels)):
        print(names[i], labels[i], cofactor_sites[i])

    # fig = plot_with_protein_from_pdb(
    # pdb_path=input_path,
    # cofactor_coords=cofactor_sites,
    # labels=labels,
    # names=names,
    # out_html="protein_plot.html"
    # )

        

    # go from active site moving away
    # check surrounding and find matching cofactor

    # return output file
    pass
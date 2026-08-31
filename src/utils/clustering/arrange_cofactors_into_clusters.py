from typing import List
from tqdm import tqdm
from src.utils.geometric.nearest_neighbor_clustering import nn_radius_clustering

def arrange_cofactors_into_clusters(
    cofactor_sites: List[tuple[float, float, float]], nn_clustering_radius : float
) -> dict:
    """ """
    # perform clustering
    labels = nn_radius_clustering(cofactor_sites, radius=nn_clustering_radius)

    # sort into clusters
    cluster = {}
    for i in tqdm(range(len(labels))):
        if labels[i] not in cluster.keys():
            cluster[labels[i]] = []
        cluster[labels[i]].append(cofactor_sites[i])

    return cluster
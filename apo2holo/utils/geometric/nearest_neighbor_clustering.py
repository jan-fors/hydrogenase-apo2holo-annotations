from sklearn.neighbors import NearestNeighbors
from scipy.sparse.csgraph import connected_components

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
"""
"""
# imports
from collections import Counter
import networkx as nx
import numpy as np
from Bio.PDB import NeighborSearch
from sklearn.feature_extraction import FeatureHasher
from pathlib import Path
from apo2holo.pdb.pdb_handler import get_structure

def create_wl_hash(structure_path: Path, point: tuple, fingerprint_radius: float, params) -> np.array:
    """
    """
    point = np.asarray(point, dtype=float)
    n_wl_columns = params["n_wl_columns"]
    wl_iters = params["wl_iters"]


    # create feature_hasher
    hasher = FeatureHasher(n_features=n_wl_columns, input_type="dict", alternate_sign=False)

    # get structure model
    model = get_structure(structure_path)

    # make neighborsearch
    ns = _make_neighbor_search(model)

    # create fingerprint
    fingerprint = sphere_wl_columns(hasher, ns, point, fingerprint_radius, wl_iters)

    return fingerprint

def _sphere_atom_graph(ns, center, r=10.0, contact_cut=4.5, bond_cut=2.1,
                      label="name", whole_residues=True):
    """
    Nodes: heavy atoms, label 'element' (C, N, O, S) or 'name' (e.g. 'ALA:CB').
    Edges: distance < bond_cut -> 'bond', < contact_cut -> 'contact' (attr 'etype').
    whole_residues=True keeps residues intact; False cuts strictly at r.
    """
    if whole_residues:
        atoms = [a for res in _sphere_residues(ns, center, r)
                 for a in res if a.element != "H"]
    else:
        hits = ns.search(np.asarray(center, dtype=float), r)
        atoms = [a for a in hits if a.get_parent().id[0] == " "]
 
    idx = {id(a): i for i, a in enumerate(atoms)}
    G = nx.Graph()
    for a in atoms:
        if label == "element":
            lab = a.element
        else:
            lab = f"{a.get_parent().get_resname()}:{a.get_id()}"
        G.add_node(idx[id(a)], label=lab)
 
    if len(atoms) > 1:
        for a, b in NeighborSearch(atoms).search_all(contact_cut):
            d = a - b  # Biopython: distance in Angstrom
            G.add_edge(idx[id(a)], idx[id(b)],
                       etype="bond" if d < bond_cut else "contact")
    return G

def _make_neighbor_search(model):
    """
    model: a Bio.PDB Model (e.g. structure[0]).
    """
    atoms = [a for a in model.get_atoms() if a.element != "H"]
    return NeighborSearch(atoms)

def _sphere_residues(ns, center, r):
    """
    Standard residues with at least one heavy atom within r of center.
    """
    hits = ns.search(np.asarray(center, dtype=float), r)
    residues = {id(a.get_parent()): a.get_parent() for a in hits}.values()
    return [res for res in residues if res.id[0] == " "]  # drop water/ligands

def wl_features(G, node_attr, edge_attr, iters=2):
    """
    Counts of WL patterns, tagged with their iteration (0 = node labels).
    """
    hashes = nx.weisfeiler_lehman_subgraph_hashes(
        G, node_attr=node_attr, edge_attr=edge_attr, iterations=iters)
    c = Counter(f"0_{G.nodes[n][node_attr]}" for n in G)
    for per_node in hashes.values():
        c.update(f"{i + 1}_{h}" for i, h in enumerate(per_node))
    return c
 
 
def wl_vector(hasher, G, node_attr, edge_attr, iters=2):
    """
    Fixed-length vector of K integer counts.
    """
    x = hasher.transform([wl_features(G, node_attr, edge_attr, iters)])
    return x.toarray()[0].astype(int)

def sphere_wl_columns(hasher, ns, center, r=10.0, iters=2, prefix="wl_"):
    """
    Returns {'wl_0': int, ..., 'wl_{K-1}': int} for one sphere.
    """
    G = _sphere_atom_graph(ns, center, r=r)
    print(G)
    
    vec = wl_vector(hasher, G, node_attr="label", edge_attr="etype", iters=iters)

    return np.array(vec)#{f"{prefix}{i}": int(v) for i, v in enumerate(vec)}

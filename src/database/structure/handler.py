from pathlib import Path
from src.io.writers.printl import printl
from src.database.structure.StructureDB import StructureDB
from src.pdb.pdb_handler import extract_protein_chains_from_file, extract_cofactors
from tqdm import tqdm
from src.parser.parse_msearch_output import parse_msearch_output
from src.filter.filter_msearch_output import filter_msearch_output
from typing import List
from src.parser.parse_t import parse_t
from src.parser.parse_u import parse_u
import os
from src.filter.apply_whitelist import apply_active_site_whitelist, apply_fes_whitelist
from src.utils.geometric.apply_transformations_to_cofactors import (
    apply_transformations_to_cofactors,
)
from src.utils.geometric.calculate_geometric_centers import calculate_geometric_centers
from src.utils.clustering.arrange_cofactors_into_clusters import (
    arrange_cofactors_into_clusters,
)


def identify_possible_active_site_pockets(
    structure_db_path: Path,
    input_path: Path,
    out: Path,
    chain_dir: Path,
    nn_clustering_radius: float,
    fident_threshold: float,
    bits_threshold: float,
) -> tuple[list, dict]:
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
    chain_paths = extract_protein_chains_from_file(input_path=input_path, out=out)

    # search against active-site structure db
    printl("Search against active-site structure database for structural homologs ...")
    active_site_sites = _active_site_structure_db_search(
        structureDB=active_site_structureDB,
        chain_paths=chain_paths,
        out=out,
        chain_dir=chain_dir,
        fident_threshold=fident_threshold,
        bits_threshold=bits_threshold
    )

    # intermediate check
    if len(active_site_sites) == 0:
        print("No possible active sites found")
        return None, None
    else:
        printl(f"Found {len(active_site_sites)} possible active-site pockets")

    # sort into clusters
    printl("Cluster cofactor sites into possible binding pockets ...")
    clusters = arrange_cofactors_into_clusters(cofactor_sites=active_site_sites, nn_clustering_radius=nn_clustering_radius)
    printl(f"Sorted possible pockets into {len(clusters)} clusters.")

    return chain_paths, clusters


def _active_site_structure_db_search(
    structureDB: StructureDB, chain_paths: List[Path], out: Path, chain_dir: Path, fident_threshold :float, bits_threshold : float
) -> List[tuple[float, float, float]]:
    """ """

    # perform foldseek searches
    fd_res = []
    for cp in tqdm(chain_paths):
        fd_res.append(structureDB.search(cp, out))

    # parse the results
    active_site_sites = []
    names = []  # ONLY relevant if plot is created

    printl(
        "Parse the results of the active site homology search in order to identify possible binding pockets"
    )
    for fd_r in tqdm(fd_res):
        # parse output
        r = parse_msearch_output(fd_r)

        # apply filters
        subset = filter_msearch_output(r, fident_threshold, bits_threshold)

        # extract coordinates for cofactors
        for i in range(subset.shape[0]):
            target = subset.iloc[i]["target"]
            u = subset.iloc[i]["u"]
            u_vec = parse_u(u)

            t = subset.iloc[i]["t"]
            t_vec = parse_t(t)

            # get structure path
            structure_path = structureDB.get_structure_path(target, chain_dir)
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


def _fes_structure_db_search(
    structureDB: StructureDB, chain_paths: List[Path], out: Path, chain_dir: Path, fident_threshold, bits_threshold
) -> List[tuple[float, float, float]]:
    """ """
    # perform foldseek searches
    fd_res = []
    for cp in tqdm(chain_paths):
        fd_res.append(structureDB.search(cp, out))

    # parse the results
    fes_sites = []
    names = []  # ONLY relevant if plot is created

    printl(
        "Parse the results of the fes homology search in order to identify possible binding pockets"
    )
    for fd_r in tqdm(fd_res):
        # parse output
        r = parse_msearch_output(fd_r)

        # apply filters
        subset = filter_msearch_output(r, fident_threshold, bits_threshold)

        printl(f"{subset.shape[0]} foldseek hits after filtering")

        # extract coordinates for cofactors
        for i in range(subset.shape[0]):
            target = subset.iloc[i]["target"]
            u = subset.iloc[i]["u"]
            u_vec = parse_u(u)

            t = subset.iloc[i]["t"]
            t_vec = parse_t(t)

            # get structure path
            structure_path = structureDB.get_structure_path(target, chain_dir)
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


def identify_possible_fes_pockets(
    structure_db_path: Path,
    out: Path,
    chain_dir: Path,
    chain_paths,
    nn_clustering_radius: float,
    fident_threshold: float,
    bits_threshold: float,
):
    """ """
    # initialize database
    printl("Initialize fes structure database ...")
    fes_structureDB = StructureDB()
    fes_structure_db_path = structure_db_path
    fes_structureDB.load(structure_db_path=fes_structure_db_path)

    # search against structure db
    printl("Search against fes structure database for structural homologs ...")
    fes_sites = _fes_structure_db_search(
        structureDB=fes_structureDB,
        chain_paths=chain_paths,
        out=out,
        chain_dir=chain_dir,
        fident_threshold=fident_threshold,
        bits_threshold=bits_threshold
    )

    # intermediate check
    if len(fes_sites) == 0:
        print("No possible fes found")
        return None
    else:
        printl(f"Found {len(fes_sites)} possible fes pockets")

    # sort into clusters
    printl("Cluster cofactor sites into possible binding pockets ...")
    clusters = arrange_cofactors_into_clusters(cofactor_sites=fes_sites, nn_clustering_radius=nn_clustering_radius)
    printl(f"Sorted possible pockets into {len(clusters)} clusters.")

    return clusters

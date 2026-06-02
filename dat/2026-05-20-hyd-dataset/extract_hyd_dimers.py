"""
Script that takes a pdb file containing a nife hydrogenase dimer
and extracts it into a separate pdb file containing only the dimer. 
The script also renames the file to rcsb_chains.pdb.

Differentiates between two cases
1) large subunit and small subunit hydrogenases
2) alpha beta gamma subunit hydrogenases, like most F420 types
"""

import argparse
import os
from Bio import SeqIO
import pyhmmer
from pyhmmer.plan7 import HMMFile
from pyhmmer import easel, plan7
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Structure
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress only Biopython PDB construction warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

def _classify_sequence(sequence, hmm_list):
    """
    Classify a sequence by scoring it against a list of HMMs.
    Returns the label of the best hit, or None if no hits.
    """
    # Clean the sequence text
    clean_seq = str(sequence).strip().upper().replace("*", "")
    if not clean_seq:
        return None

    alphabet = easel.Alphabet.amino()

    # Create the DigitalSequence properly
    text_seq = easel.TextSequence(
        name=b"my_query_sequence", 
        sequence=clean_seq
    )
    query_seq = text_seq.digitize(alphabet)

    # 1. CRITICAL FIX: Wrap your sequence inside a DigitalSequenceBlock
    # This matches the explicit C-level signature that Cython expects.
    sequence_block = easel.DigitalSequenceBlock(alphabet, [query_seq])

    pipeline = plan7.Pipeline(alphabet=alphabet, E=1e-5)
            
    best_hit = None
    best_score = float('-inf')
    
    for hmm_label, hmm_model in hmm_list:
        # 2. Pass the sequence_block here instead of a raw python list []
        hits = pipeline.search_hmm(hmm_model, sequence_block)
        
        if len(hits) > 0:
            top_hit_score = hits[0].score 
            
            if top_hit_score > best_score:
                best_score = top_hit_score
                best_hit = hmm_label
    
    return best_hit

def _calculate_chain_mass_center(pdb_file, chain_id):
    parser = PDBParser()
    # ... your existing PDBParser setup ...
    structure = parser.get_structure("id", pdb_file)
    
    # FIX: Don't hardcode [0]. Grab the first available model dynamically, 
    # whether its ID is 0, 1, or a string.
    try:
        model = next(structure.get_models())
    except StopIteration:
        raise ValueError(f"The PDB file {pdb_file} contains no models.")

    # Check if the chain actually exists in this model before grabbing it
    if chain_id not in model:
        # Debugging aid: Show what chains ARE actually available to help you troubleshoot
        available_chains = [c.id for c in model]
        raise KeyError(
            f"Chain '{chain_id}' not found in the first model of {pdb_file}. "
            f"Available chains are: {available_chains}"
        )
        
    chain_object = model[chain_id]
    coord_sum = [0.0, 0.0, 0.0]
    atom_count = 0
    
    # Calculate mass center directly using the atoms already in memory
    for atom in chain_object.get_atoms():
        # You can optionally weight by element mass here, 
        # or do a simple geometric center:
        coord_sum[0] += atom.coord[0]
        coord_sum[1] += atom.coord[1]
        coord_sum[2] += atom.coord[2]
        atom_count += 1
        
    if atom_count == 0:
        raise ValueError(f"Chain {chain_object.id} has no atoms!")
        
    return [c / atom_count for c in coord_sum]

def _extract_submer(pdb_file, chains, output_path):
    """
    Extract a submer from a pdb file containing the specified chains and save it to an output file.
    """
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    submer_structure = structure.copy()

    for model in submer_structure:
        all_chain_ids = [chain.id for chain in model]
        
        for chain_id in all_chain_ids:
            if chain_id not in chains:
                model.detach_child(chain_id)

        single_model_structure = Structure.Structure("single_model")
        single_model_structure.add(model.copy())
        break
    
    io = PDBIO()
    io.set_structure(single_model_structure)
    name = pdb_file.split("/")[-1].split(".")[0]
    conc_chains = "".join(chains)
    output_file = os.path.join(output_path, f"{name}_{conc_chains}.pdb")
    io.save(output_file)

def extract_dimer(pdb_file, output_dir):
    """
    Read all chains from the input file and extract the sequences
    score all sequences using hmms and sort into small and large subunits, or alpha beta gamma
    match subunits next to each other and extract the dimer
    save dimers in output dir
    """
    # read all chains from the input file
    chain_sequence_pairs = []

    for record in SeqIO.parse(pdb_file, "pdb-seqres"):
        print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        chain_sequence_pairs.append((record.annotations["chain"], record.seq))

    if len(chain_sequence_pairs) == 0:
        print("No chains found in the input file.")
        return
    
    if len(chain_sequence_pairs) == 1:
        print("Only one chain found in the input file. No dimer to extract.")
        return
    
    # classify sequences using hmms
    with HMMFile("hmms/NF007550.0.hmm") as hmm_file:
        hmm_large = hmm_file.read()

    with HMMFile("hmms/NF045520.1.hmm") as hmm_file:
        hmm_small = hmm_file.read()

    with HMMFile("hmms/TIGR03295.1.hmm") as hmm_file:
        hmm_alpha = hmm_file.read()

    with HMMFile("hmms/TIGR03289.1.hmm") as hmm_file:
        hmm_beta = hmm_file.read()

    with HMMFile("hmms/TIGR03294.1.hmm") as hmm_file:
        hmm_gamma = hmm_file.read()


    hmm_list = [
        ("large_subunit", hmm_large), 
        ("small_subunit", hmm_small), 
        ("alpha_subunit", hmm_alpha), 
        ("beta_subunit", hmm_beta), 
        ("gamma_subunit", hmm_gamma)]
    
    # Get the chains that actually have coordinates in the structure
    parser_check = PDBParser()
    structure_check = parser_check.get_structure("id", pdb_file)
    model_check = next(structure_check.get_models())
    available_chains = {chain.id for chain in model_check}

    classified_sequences = []
    for chain, sequence in chain_sequence_pairs:
        if chain not in available_chains:
            print(f"Skipping chain {chain}: present in SEQRES but has no ATOM records.")
            continue
        print(f"Classifying chain {chain} with sequence {sequence}")
        label = _classify_sequence(str(sequence), hmm_list)
        classified_sequences.append((chain, sequence, label))
        print(f"Chain {chain}: {label}")

    # calculate chain mass centers
    with_mass_centers = []
    for chain, sequence, label in classified_sequences:
        print(f"Calculating mass center for chain {chain} with label {label}")
        mass_center = _calculate_chain_mass_center(pdb_file, chain)
        
        print(f"Mass center for chain {chain}: {mass_center}")
        with_mass_centers.append((chain, sequence, label, mass_center))

    # sort and export dimers/trimers if f420
    large_subunits = [item for item in with_mass_centers if item[2] == "large_subunit"]
    small_subunits = [item for item in with_mass_centers if item[2] == "small_subunit"]
    alpha_subunits = [item for item in with_mass_centers if item[2] == "alpha_subunit"]
    beta_subunits = [item for item in with_mass_centers if item[2] == "beta_subunit"]
    gamma_subunits = [item for item in with_mass_centers if item[2] == "gamma_subunit"]

    for large in large_subunits:
        smallest_distance = float('inf')
        closest_small = None
        for small in small_subunits:
            # check if mass centers are close enough to be a dimer
            distance = ((large[3][0] - small[3][0]) ** 2 + (large[3][1] - small[3][1]) ** 2 + (large[3][2] - small[3][2]) ** 2) ** 0.5
            print(f"Distance between large {large[0]} and small {small[0]}: {distance}")
            
            if distance < smallest_distance:
                smallest_distance = distance
                closest_small = small

        if closest_small is not None:
            print(f"Closest small subunit to large subunit {large[0]} is {closest_small[0]} with distance {smallest_distance}")
            # extract dimer and save to output dir
            
            _extract_submer(pdb_file, [large[0], closest_small[0]], output_dir)


    for alpha in alpha_subunits:
        closest_beta = None
        closest_gamma = None
        smallest_beta_distance = float('inf')
        smallest_gamma_distance = float('inf')

        for beta in beta_subunits:
            distance = ((alpha[3][0] - beta[3][0]) ** 2 + (alpha[3][1] - beta[3][1]) ** 2 + (alpha[3][2] - beta[3][2]) ** 2) ** 0.5
            print(f"Distance between alpha {alpha[0]} and beta {beta[0]}: {distance}")
            
            if distance < smallest_beta_distance:
                smallest_beta_distance = distance
                closest_beta = beta

        for gamma in gamma_subunits:
            distance = ((closest_beta[3][0] - gamma[3][0]) ** 2 + (closest_beta[3][1] - gamma[3][1]) ** 2 + (closest_beta[3][2] - gamma[3][2]) ** 2) ** 0.5
            print(f"Distance between beta {closest_beta[0]} and gamma {gamma[0]}: {distance}")
            
            if distance < smallest_gamma_distance:
                smallest_gamma_distance = distance
                closest_gamma = gamma

        if closest_beta is not None and closest_gamma is not None:
            print(f"Closest beta subunit to alpha subunit {alpha[0]} is {closest_beta[0]} with distance {smallest_beta_distance}")
            print(f"Closest gamma subunit to beta subunit {closest_beta[0]} is {closest_gamma[0]} with distance {smallest_gamma_distance}")
            # extract trimer and save to output dir TODO
            _extract_submer(pdb_file, [alpha[0], closest_beta[0], closest_gamma[0]], output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract nife hydrogenase dimers from pdb files.")
    parser.add_argument("pdb_file", type=str, help="Path to the input pdb file containing the dimer.")
    parser.add_argument("output_dir", type=str, help="Directory where the extracted dimer will be saved.")
    
    args = parser.parse_args()
    
    extract_dimer(args.pdb_file, args.output_dir)
"""
Takes the folder 2026-03-06_dedup-dimers as input and 
extracts the sequences to a fasta file. The output is written to 2026-03-06_dedup-dimers.fasta.
Both chains are extracted and sorted by length, then the sequence is written like this:
<longer_chain_sequence>:<shorter_chain_sequence>
"""

import os
from Bio import SeqIO
input_folder = "2026-03-06-dedup_dimers"
output_file = "2026-03-06_dedup-dimers.fasta"

with open(output_file, "w") as out_f:
    for filename in os.listdir(input_folder):
        if filename.endswith(".pdb"):
            filepath = os.path.join(input_folder, filename)
            # Extract sequences from the PDB file
            sequences = {}
            for record in SeqIO.parse(filepath, "pdb-seqres"):
                sequences[record.id] = str(record.seq)
            # Sort chains by length
            sorted_chains = sorted(sequences.items(), key=lambda x: len(x[1]), reverse=True)
            if len(sorted_chains) >= 2:
                longer_chain_seq = sorted_chains[0][1]
                shorter_chain_seq = sorted_chains[1][1]
                out_f.write(f">{filename}\n{longer_chain_seq}:{shorter_chain_seq}\n")
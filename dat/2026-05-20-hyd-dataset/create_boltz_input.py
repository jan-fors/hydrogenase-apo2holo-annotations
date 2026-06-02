"""

"""
import os
import argparse
from pathlib import Path
from Bio import SeqIO 

LETTERS = ["A", "B", "C", "D"]

def create_boltz_input(sequences_fasta : Path, msa_dir : Path, output_dir : Path):
    """"""
    # read all files from the msa dir
    msa_dir_files = os.listdir(msa_dir)
    msa_dir_files = [int(x.split(".")[0]) for x in msa_dir_files]
    msa_dir_files = sorted(msa_dir_files)
    

    counter = 0
    hyds = {}
    # read all sequences
    for seq_record in SeqIO.parse(sequences_fasta, "fasta"):
        seq_record_id = seq_record.id
        rcsb = seq_record_id.split("_")[0]
        
        if rcsb not in hyds.keys():
            hyds[rcsb] = []

        hyds[rcsb].append({
            "sequence": str(seq_record.seq),
            "msa": Path("/home/jan/dat/2026-05-20-hyd-dataset/unique_structures_by_tm/msas") / Path(str(msa_dir_files[counter]) + ".a3m")
            })
        
        counter += 1

    for h in hyds.keys():
        res = "sequences:\n"

        for i in hyds[h]:
            res += "  - protein:\n"
            res += f"      id: [{LETTERS[hyds[h].index(i)]}]\n"
            res += f"      sequence: {i['sequence']}\n"
            res += f"      msa: {i['msa']}\n"

        with open(output_dir/ Path(h + "_boltz_input.yaml"), "w") as f:
            f.write(res)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("sequences", type=Path, help="Path to the sequences fasta")
    parser.add_argument("msa_dir", type=Path, help="Path to the directory containing msas")
    parser.add_argument("output", type=Path, help="Path to the output directory")

    args = parser.parse_args()

    if not os.path.exists(args.sequences):
        raise ValueError("Sequences file does not exist")
    
    if not os.path.exists(args.msa_dir):
        raise ValueError("MSA dir does not exist")
    
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    create_boltz_input(args.sequences, args.msa_dir, args.output)
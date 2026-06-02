import sys
import pandas as pd

def count_representatives(run_name: str) -> int:
    tsv_path = f"{run_name}_cluster.tsv"
    df = pd.read_csv(tsv_path, sep="\t", header=None, names=["representative", "member"])
    n_reps = df["representative"].nunique()
    return n_reps

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_reps.py <run_name>")
        sys.exit(1)

    run_name = sys.argv[1]
    n = count_representatives(run_name)
    print(f"{n} representatives in {run_name}_cluster.tsv")
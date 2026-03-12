import os
import random
import shutil
from pathlib import Path

source_dir = Path("dat/dedup_dimer")
target_dir = Path("dat/benchmarking_subsets")
n_sets = 10

# alle Dateien sammeln
files = [f for f in source_dir.iterdir() if f.is_file()]

# zufällig mischen
random.shuffle(files)

# Unterordner erstellen
set_dirs = []
for i in range(n_sets):
    d = target_dir / f"set_{i+1}"
    d.mkdir(parents=True, exist_ok=True)
    set_dirs.append(d)

# Dateien möglichst gleichmäßig verteilen
for i, file in enumerate(files):
    dest = set_dirs[i % n_sets] / file.name
    shutil.copy2(file, dest)

print("Fertig.")
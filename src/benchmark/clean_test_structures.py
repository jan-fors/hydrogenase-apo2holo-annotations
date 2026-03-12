

from pathlib import Path
import os

folder = Path("/home/solar/Documents/hiwi/deephyds/hydrogenase-apo2holo-annotations/db/benchmarking_dbs/set_2/test")

for pdb_file in folder.glob("*.pdb*"):
    temp_file = pdb_file.with_suffix(".tmp")

    with open(pdb_file, "r") as fin, open(temp_file, "w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                fout.write(line)
            elif line.startswith(("TER", "END")):
                fout.write(line)

    os.replace(temp_file, pdb_file)

print("Fertig.")
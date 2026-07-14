from pathlib import Path

folder = Path("/home/jan/hydrogenase-apo2holo-annotations/dat/2026-05-26-cof-dataset/cofactors_of_interest/F3S/pdb")

for pdb_file in folder.glob("*.pdb"):
    with open(pdb_file) as f:
        lines = f.readlines()

    output = []
    in_first_model = False
    has_models = False

    for line in lines:
        if line.startswith("MODEL"):
            has_models = True
            if not in_first_model:
                in_first_model = True
                output.append(line)
            continue

        if line.startswith("ENDMDL"):
            output.append(line)
            break

        if in_first_model:
            output.append(line)

    # If there were no MODEL records, keep the file unchanged
    if not has_models:
        continue

    with open(pdb_file, "w") as f:
        f.writelines(output)

print("Done.")
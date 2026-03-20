# hydrogenase-apo2holo-annotations

This repository provides a tool to annotate apo **nife** hydrogenase structures with
information required to obtain the corresponding holoenzyme state.

## Scope
- Identify required cofactors (e.g. metal clusters)
- Assign cofactor types to binding sites
- Provide residue- and position-level placement information

## Input
- PDB or mmCIF file containing the hydrogenase apo-enzyme

## Output
- Structured annotation data (JSON / tabular)
- No structure generation or modification
- Input for specific structure prediction models (e.g. yaml for boltz-2)


## Run
### apo2holo
```sh
python -m src.cli <INPUT STRUCTURE> <RUN_NAME> -o <OUTPUT_DIR>
```
### Build 
```sh
python -m src.build <RAW_STRUCTURE_DIR> <OUTPUT_DIR FOR CHAINS> --structure-db-path <PATH> --fingerprint-db-path <PATH>
```


# References
- Rodrigues JPGLM, Teixeira JMC, Trellet M and Bonvin AMJJ.
pdb-tools: a swiss army knife for molecular structures. 
F1000Research 2018, 7:1961 (https://doi.org/10.12688/f1000research.17456.1) 
- https://github.com/steineggerlab/foldseek?tab=readme-ov-file#search
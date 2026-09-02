# hydrogenase-apo2holo-annotations

This repository provides a tool to annotate apo **nife** hydrogenase structures with
information required to obtain the corresponding holoenzyme state.

## Scope
- Identify required cofactors (e.g. metal clusters)
- Assign cofactor types to binding sites
- Provide residue- and position-level placement information

## Input
- PDB file containing the hydrogenase apo-enzyme

## Output
- Structured annotation data (JSON / tabular)
- No structure generation or modification
- Input for specific structure prediction models (e.g. yaml for boltz-2)

## Run
In order to use the program create a conda env with:
```sh
conda env create -f environment.yml
```
### apo2holo
Before running the script for the first time check `src/utils/constants.py`. 
In order for the program to run correctly, the paths to the databases have to be set correctly.
```sh
python -m src.cli <INPUT STRUCTURE> <RUN_NAME> -o <OUTPUT_DIR>
```

*Example*: `python -m src.cli example/apo_3RGW.pdb -o out testi --structure-db-path db/structureDB/structureDB --fingerprint-db-path db/fingerprintDB`
### Build Databases
```sh
python -m src.build <RAW_STRUCTURE_DIR> <OUTPUT_DIR FOR CHAINS> --structure-db-path <PATH> --fingerprint-db-path <PATH>
```
- `OUTPUT_DIR FOR CHAINS`: Each structure is split into its chains and the chains are copied into this folder,
- `--structure-db-path`: Path for the final structure database. Database is created inside a new folder inside of this specified path 

*Example*: `python -m src.build dat/2026-03-06-dedup_dimers db/single_chains --structure-db-path db/structureDB --fingerprint-db-path db/fingerprintDB`

# References
- Rodrigues JPGLM, Teixeira JMC, Trellet M and Bonvin AMJJ.
pdb-tools: a swiss army knife for molecular structures. 
F1000Research 2018, 7:1961 (https://doi.org/10.12688/f1000research.17456.1) 
- https://github.com/steineggerlab/foldseek?tab=readme-ov-file#search   


#pip install git+https://github.com/you/your-repo.git
# hydrogenase-apo2holo-annotations

This repository provides tools to annotate apo hydrogenase structures with
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
### Build structure DB
```sh
python -m src.db.build_structure_db <DIR WITH PDB STRUCTURES> <DIR FOR CLEANED STRUCTURE CHAINS>
```
- creates structureDB



# References
- Rodrigues JPGLM, Teixeira JMC, Trellet M and Bonvin AMJJ.
pdb-tools: a swiss army knife for molecular structures. 
F1000Research 2018, 7:1961 (https://doi.org/10.12688/f1000research.17456.1) 
- 
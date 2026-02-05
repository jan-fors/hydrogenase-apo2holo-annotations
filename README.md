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



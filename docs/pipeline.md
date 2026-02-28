# Pipeline

Follow the procedure described in AlphaFill Paper with slight differences.

1. Amino-acid sequence of structure is BLASTed against the sequence file of the database, which contains known hydrogenase structures and proteins that contain FeS-Cofactors. Return hits of database sorted by E values.
2. Structurally align the hits with the input protein on the Ca-atoms of the residues that match in the BLAST alignment.

-> 1 and 2 with foldseek

3. Identify the positions of the Active Site as well as the FeS Cofactors
4. Group the cofactors into the groups: activesite, proximal, medial and distal cluster



5. Check the surrounding of the mass center and create aminoacid fingerprints
6. select the most reasonable combination of cofactors based on the fingerprints.
7. create output file and return
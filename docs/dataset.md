# Datasets

## 2026-03-06-dataset
To aquire the structures for the initial dataset structures by search for `nife hydrogenase` on rcsb -> (548 structures 06.03.26).
In a first step the large heteromultimer complexes where split up into heterodimers.
The Set of heterodimers was then deduplicated by performing a series of step using `mmseqs` (see `dat/2026-03-06-dataset/2026-03-06-unique_hyd_structured_creation`)
The final dataset `2026-03-06-test-data` contains 27 distinct nife hydrogenases.

## 2026-05-08-dataset
Even if the application of this script is for nife hydrogenases a large amount of different structures is needed in order to accurately predict the correct cofactor for each pocket.
Therefore experimental structures for each Cofactor of interest have been downloaded and included in the dataset.
For the cofactors except 4Fe4S every structure was taken. For 4Fe4S the structures where filtered previously by: experimental, Protein, Oxidoreductase and newer than 2020-01-01.

All the structures were clustered using `foldseek easy-multimercluster` with the following params: `-c 0.8 --cov-mode 0 --min-seq-id 1`.
The resulting dataset `2025-05-08-train-data` contains 543 distinct structures.

## Intersection of both datasets
=== IN BOTH (intersection) ===
3usc.pdb
5a4f.pdb
5a4i.pdb
5adu.pdb
5jrd.pdb
5lry.pdb
6g7m.pdb
7nem.pdb
7utd.pdb
7uus.pdb
9nez.pdb
9r52.pdb
9r6z.pdb

=== IN A BUT NOT IN B ===
3ayz.pdb
4ci0.pdb
4omf.pdb
5aa5.pdb
5lmm.pdb
5odq.pdb
5xvb.pdb
5xvc.pdb
6ehq.pdb
6gam.pdb
6qgt.pdb
6qii.pdb
9erb.pdb
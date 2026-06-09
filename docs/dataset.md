# Datasets

## Hydrogenase Dataset
### Experimental Data
Constructed by appying the following pipeline:
1. Search rcsb.org for **hydrogenases** with filters set on "experimental", "protein", and "oxidoreductase". (20.05.26 693 hits)
2. Batch download
3. Extracting NiFe-Hydrogenase into `initial_structures`
    1. `unzip '*.zip'`
    2. `gunzip *.gz`
4. Extract corresponding Hydrogenase Dimers using `extract_hyd_dimers.py` into `functional_units`
5. Deduplicating dataset by a tm-score of > 95% using `foldseek easy-multimercluster functional_units/ foldseek/clu_tm095 tmp --multimer-tm-threshold 0.95 --cov-mode 0`, representatives are copied into `unique_structures_by_tm`

```dataset_summary.py
Structures 10
========== FeS Cluster ==========
SF4: 27
F3S: 6
SF3: 1
```

6. Expand dataset by clustering the cofactor fingerprints inside of each cluster with `expand_by_fingerprint_clustering.py`, Faktor=1.0 -> saved into `unique_structures_after_cl_F1`

```dataset_summary.py
Structures 38
========== FeS Cluster ==========
SF4: 74
F3S: 32
F4S: 6
SF3: 5
```

### Apoenzyme Structures
1. The sequences of the `unique_structures_after_cl_F1`/`unique_structures_by_tm` were downloaded from rcsb.
2. MSAs of the single chains were made using `colabfold_search` against `uniref30_2302`
3. Boltz-2 Input yamls were generated using `create_boltz_input.py <sequences fasta> <msa dir> <output>`

## Cofactor Dataset
cofactors of interest: F3S, FSX, SF4, F4S, SF3, ER2
1. Search rcsb for all exp. structures that contain any of those cofactors of interest and download. Filter: `protein`
    - F3S (372 hits)
    - FSX (4 hits) -> *use?*
    - SF4 (2505 hits -> 2113 after download)
    - F4S (15 hits)
    - SF3 (33 hits)
    - ER2 (3 hits) -> *use?*
2. deduplicate all structures by using `foldseek easy-multimersearch F3S/pdb F4S/pdb SF3/pdb SF4/pdb foldseek/<result> foldseek/tmp --min-seq-id 1.00 --cov-mode 0` leaving 742 stuctures

```dataset_summary.py
Structures 1329
========== FeS Cluster ==========
SF4: 4162
F3S: 196
SF3: 22
F4S: 8
```


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

SF4: 1349
F3S: 171
SF3: 8
F4S: 3
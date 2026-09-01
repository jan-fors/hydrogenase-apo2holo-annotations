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

6. Expand dataset by clustering the cofactor fingerprints (Aminoacid Count) inside of each cluster with `expand_by_fingerprint_clustering.py`, Faktor=1.0 -> saved into `unique_structures_after_cl_F1`

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
cofactors of interest: F3S, SF4, F4S, SF3
1. Search rcsb for all exp. structures that contain any of those cofactors of interest and download. Filter: `protein`
    - F3S (372)
    - SF4 (753), additional filter were applied: {asymmetric unit, oxidoreductase, experimental}
    - F4S (15)
    - SF3 (33)

    Copy all of them into the folder `structures`

2. deduplicate all structures by using `foldseek easy-multimercluster structures/ foldseek/cof_dataset foldseek/tmp --min-seq-id 1.00 --cov-mode 0` leaving 584 stuctures
3. These are copied into `seq_repr_v3` and used as training data

```dataset_summary.py
Structures 584
========== FeS Cluster ==========
SF4: 1835
F3S: 308
SF3: 24
F4S: 9
```
# TODOs
## General
- [X] finalize output
- [X] boltz output
- [ ] change raw output so that they are linear distanced from active site
- [ ] check whether smiles are useful for boltz
- [X] remove pdb files etc after running the program

## Benchmarking
- [X] update benchmarking engine
    - [X] given a *complete benchmark cv* folder -> allow to:
        - [X] train different models
        - [X] run program against different models
        - [X] compare runs
- [ ] analyze results and find ways to increase accuracy
    - [ ] different dataset
        - [ ] TEST: oversample 3Fe4S cluster or undersample 4Fe4S
            1. deduplicate by fingeprints
            2. oversample using `random_oversampler.py`
        - [ ] TEST: generally less items
    - [ ] different training
    - [ ] different models for FeS cluster and acitve site
    - [ ] ignore all unimportant aminoacids? -> just CYS, HIS, etc.. 
    - [ ] when creating the training set -> smaller min-seq-id threshold

- [X] run best models on larger dataset
    === TOP ===
modelname f_radii search_type model_type db_R  model_number  accuracy  coverage
A     6.5          mc        rfc R6.0             3  0.933333       1.0
B     6.0          mc        mlp R6.0            19  0.933333       1.0
C     6.5          mc        rfc R6.0            11  0.933333       1.0
D     6.5         sum        lrc R6.0             2  0.933333       1.0
E     6.5          mc        rfc R6.0            10  0.933333       1.0

F     6.5         sum        mlp R6.0             8  0.900000       1.0
G     6.0         sum        mlp R6.0            12  0.900000       1.0
H     6.5          mc        svm R6.5             9  0.900000       1.0

I     6.0          mc        rfc R6.0             7       0.9       1.0        1.0
J     6.5         sum        mlp R6.0             8       0.9       0.9        1.0
K     5.0          mc        mlp R5.0             7       0.8       0.8        1.0

## IMPORTANT
- [ ] **when** are **which cofactors** removed or flagged?
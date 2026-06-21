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
        - [ ] TEST: generally less items *overfit*?
    - [ ] different training
    - [ ] different models for FeS cluster and acitve site
    - [X] ignore all unimportant aminoacids? -> just CYS, HIS, etc.. -> *no positive effect* 
    - [X] when creating the training set -> smaller min-seq-id threshold

### IMPORTANT
- [ ] test with apostructures from experimental files
- [X] remove GLY from search? -> *no positive effect*
- [ ] increase background class size? -> make background class stronger
    - [ ] background class too different?
    - [X] create artificial background class samples?

- [ ] implement absolut search

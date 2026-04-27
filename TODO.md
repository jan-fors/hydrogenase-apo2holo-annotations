# TODOs
## General
- [X] database
    - [X] fingerprint db
    - [X] structure db -> cleaned structures
- [X] broad structure
- [X] foldseek impl
- [X] identify positions
- [X] create fingerprints
- [X] group
- [ ] select
    - [X] Logistic Regression for each cluster? -> maybe later first just with similarity
- [X] split models in build db
- [X] Create FingerprintDB class
    - [X] fill
- [X] return
- [X] load fingerprintDB once and use not every time search against is called
- [X] Benchmarking überarbeiten
- [ ] finalize output
- [ ] implement models
    - [ ] svm mc & sum
    - [ ] mlp classifier mc & sum
    - [ ] pytorchnn mc & sum
- [ ] perform benchmarking grid search
    - [ ] define parameters

## Extra
- [ ] option to return top x different with scores

## Benchmarking
- [ ] colabfold structures as input for benchmarking? Loss ap02holo(colabfold prediction) against experimental prediction

## Issues:
- [X] some files have a different name inside of the structureDB therefore they are skipped -> make everyone used
- [X] look at wrong results of clustering
- [ ] too many with non conform fes cluster amount
    - [ ] bad pre filtering?
    - [X] correct subunits in teststructures? -> seems so
    - [X] active site identification is bad
    - [ ] increase fingerprint radius for distal cluster?
    - [X] **one proba search per mass center for each cluster spot?** -> not so good
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
- [X] select
    - [X] Logistic Regression for each cluster? -> maybe later first just with similarity
- [X] split models in build db
- [X] Create FingerprintDB class
    - [X] fill
- [X] return
- [X] load fingerprintDB once and use not every time search against is called
- [X] Benchmarking überarbeiten
- [ ] finalize output
- [X] implement models
    - [X] svm mc & sum
    - [X] mlp classifier mc & sum
- [X] refactor 
- [ ] perform benchmarking grid search
    - [ ] define parameters
- [ ] boltz output
- [ ] visualizing output
- [ ] check whether smiles are useful for boltz

## Extra
- [ ] option to return top x different with scores
- [ ] no reason to load all models at once - except ensemble prediction?

## Benchmarking
- [ ] colabfold structures as input for benchmarking? Loss ap02holo(colabfold prediction) against experimental prediction

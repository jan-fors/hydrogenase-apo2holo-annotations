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
- [X] visualizing output
    - [X] visualize the protein in a 3d interactive plot 
    - [X] show pockets -> shows all predicted pockets colored by cluster
    - [X] for each cluster show the fingerprint radius and create a little table with each cluster, then the probabilites for each cofactor and the winner
- [ ] boltz output
- [ ] check whether smiles are useful for boltz

## Benchmarking
- [ ] perform benchmarking grid search
    - [ ] define parameters
- [ ] colabfold structures as input for benchmarking? Loss ap02holo(colabfold prediction) against experimental prediction
- [ ] create testsets

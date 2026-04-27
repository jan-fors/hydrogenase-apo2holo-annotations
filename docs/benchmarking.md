# Benchmarking

The Benchmarking is described here.
A 10-fold CV is performed, comparing different types (`svm_sum`, `absolut`, ...) of cofactor prediction.

## Datasets
The original dataset `2026-03-06-dedup_dimers` is split into 10 randomized subsets.
For each of the 10 predictions one of them is reserved as `test` set and the other nine subsets are merged.
The databases and the different models are then created based on the merged subsets.
This results in 10 different databases and 10 different sets of models.

## Scoring
In order to assess the quality different scorings measured. 
Since we assume that each nife-hydrogenase has three FeS-clusters, we measure the accuracy on each of the clusters separatly. 
For each position, where a cofactor of interest is predicted, we measure:
- `Precision`
- `Recall`
- `Accuracy`
by comparing the prediction to the groundtruth with a simple 0-1 loss function (1 if correct and 0 if incorrect).
Moreover we calculate the prediction `coverage`.

## Models
All model classes which are available are tested.

## Parameter
The whole benchmarking is performed for different sets of parameters depending on the model.


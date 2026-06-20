# Benchmarking
## Model Benchmarking

## Complete Benchmarking
Procedure:
1. `create_subsets.py`
```sh
python -m src.benchmark.complete_benchmark.create_subsets dat/2026-05-20-hyd-dataset/unique_structures_by_tm/apoenzyme dat/2026-05-26-cof-dataset/seq_repr/ --k 5 --min_seq_id 0.95 --tmp cs_tmp --out benchmarking/complete_benchmark/subsets
```

2. `build_structure_dbs.py`
```sh
python -m src.benchmark.complete_benchmark.build_structure_dbs benchmarking/complete_benchmark/subsets
```

3. `build_fingerprintdb_base.py`
```sh
python -m src.benchmark.complete_benchmark.build_fingerprintdb_base benchmarking/complete_benchmark/subsets
```

3.2 `deduplicate_fingeprint_db.py`
```sh
python -m src.benchmark.complete_benchmark.deduplicate_fingeprint_db benchmarking/complete_benchmark/subsets db_R5.0.tsv
```

3.3 `select_features.py`
```sh
python -m src.benchmark.complete_benchmark.select_features benchmarking/complete_benchmark/select_features db_R5.0.tsv 
```

3.4 `resample_dbs.py`
```sh
python -m src.benchmark.complete_benchmark.resample_dbs benchmarking/complete_benchmark/subsets db_R5.0.tsv
```

4. `train_models.py`
```sh
python -m src.benchmark.complete_benchmark.train_models benchmarking/complete_benchmark/subsets db_R5.0.tsv model --model_dir_name mlp --parameter_grid --iter 20
```

5. `predict.py`
*predict using single model*
```sh
python -m src.benchmark.complete_benchmark.predict benchmarking/complete_benchmark/subsets/ --f_radius 5 --search_type sum --model mlp/model__0.pkl
```
*predicting on a set of models*
```sh
python -m src.benchmark.complete_benchmark.predict benchmarking/complete_benchmark/subsets/ --f_radius --search_type sum --model_dir lrc/R5.0
```
6. Analyze
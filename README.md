# Hydrogenase-apo2holo-Annotations
A machine learning based approach that identifies relevant cofactor pockets in NiFe-Hydrogenases and assigns fitting molecules. 

## Overview
Given a NiFe-Hydrogenase apostructure, the script initially performs a foldseek[^1] search to identify structural homologs using the dataset `seq_repr_v3` (see `docs/dataset.md` for more info) as database. The existing cofactor pockets from the structural hits are projected onto the input structure and possible candidate pockets are identified through clustering throse predicted pockets. Using a hierarchical classification process each pocket is assigned a cofactor, if it is predicted as real binding pocket. The predictions are based on inter protein fingerprints (see `docs/fingerprints.md`) and performed using a series of MLPClassifiers (see `docs/pipeline.md`). 

## Installation
Install directly from GitHub via pip in a new conda environment:
```sh
# in new environmetn
conda create -n apo2holo python=3.11 -c bioconda -c conda-forge pdb-tools foldseek
conda activate apo2holo

# install apo2holo
pip install git+https://github.com/solarflip/hydrogenase-apo2holo-annotations
```

Or clone and install from source (useful for development):
```sh
git clone https://github.com/solarflip/hydrogenase-apo2holo-annotations
cd hydrogenase-apo2holo-annotations
pip install -e .
```
## Usage
```sh
apo2holo [input] [commands]
```

### Example
```sh
apo2holo input.pdb -o out/ --plot
```

### Options
| Flag | Description | Default |
|---|---|---|
| `-o, --output_dir <path>` | Output dir path | `.` |
| `-c, --config` | Path to a custom config file | `$DATADIR/"config.yaml"`|
| `--plot` | Generate a result plot | |
| `--boltz` | Generate an output yaml for boltz[^2]||


## Citation
If you use this tool in your research, please cite:
t.b.p.

## License
[MIT](LICENSE)

[^1]: van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, and Steinegger M. Fast and accurate protein structure search with Foldseek. Nature Biotechnology, doi:10.1038/s41587-023-01773-0 (2023)

[^2]: Passaro, S., Corso, G., Wohlwend, J., Reveiz, M., Thaler, S., Somnath, V. R., Getz, N., Portnoi, T., Roy, J., Stark, H., Kwabi-Addo, D., Beaini, D., Jaakkola, T., & Barzilay, R. (2025). Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. *bioRxiv*. https://doi.org/10.1101/2025.06.14.659707
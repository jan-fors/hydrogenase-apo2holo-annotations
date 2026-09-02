from apo2holo.utils.paths import get_data_dir
import yaml
import os
from importlib.resources import files
from pathlib import Path
from apo2holo.build import build
from apo2holo.utils.fetch import fetch_file
import zipfile

def init():
    """
    """
    data_dir = get_data_dir()
    # check if everything already exists
    structure_db_path = data_dir / Path("db") /Path("structureDB")/Path("structureDB")
    model_dir = data_dir / Path("models")
    
    if os.path.exists(model_dir) and os.path.exists(structure_db_path):
        pass
    else:
        os.makedirs(model_dir, exist_ok=True)

        print(data_dir)

        # load config and fetch models
        default_config_path = files("apo2holo.config").joinpath("config.yml")
        with default_config_path.open("r") as f:
            config = yaml.safe_load(f)

        as_model_path = fetch_file("as_experiment_02.pkl", model_dir)
        config["models"]["active_site"]["path"] = str(as_model_path)

        fes_pocket_path = fetch_file("fes_pocket_experiment_02.pkl", model_dir)
        config["models"]["fes"]["pocket"]["path"] = str(fes_pocket_path)

        fes_type_path = fetch_file("fes_type_experiment_02.pkl", model_dir)
        config["models"]["fes"]["type"]["path"] = str(fes_type_path)

        # structure db
        db_dir = data_dir / "db"
        raw_structures_dir = db_dir / "raw_structures"
        raw_structures_dir.mkdir(parents=True, exist_ok=True)

        zip_path = db_dir / "structures.zip"
        marker = raw_structures_dir / ".extracted"
        if not marker.exists():
            fetch_file("structures.zip", db_dir) 
            print(f"Extracting {zip_path} to {raw_structures_dir}...")
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(raw_structures_dir)
            marker.touch()
            zip_path.unlink()  
        else:
            print("Raw structures already extracted, skipping download.")

        chain_dir = db_dir/Path("single_chains")
        os.makedirs(chain_dir, exist_ok=True)
        structure_db_path = db_dir/Path("structureDB")
        os.makedirs(structure_db_path, exist_ok=True)
        structure_db_path = structure_db_path/Path("structureDB")

        build(input_structure_dir=raw_structures_dir/Path("seq_repr_v3"),output_structure_dir=chain_dir, structure_db_path=structure_db_path,jobs=8)

        config["structure_db"]["path"] = str(structure_db_path)
        config["structure_db"]["chain_dir"] = str(chain_dir)
            
        # write config for usage
        user_config_path = data_dir / "config.yaml"
        with user_config_path.open("w") as f:
            yaml.safe_dump(config, f)
            
        print(f"Config written to {user_config_path}")
from pathlib import Path
import yaml

def read_config(config_path : Path) -> dict:
    """
    """
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    return config
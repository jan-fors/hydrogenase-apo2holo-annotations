from pathlib import Path
from platformdirs import user_data_dir

APP_NAME = "apo2holo"

def get_data_dir() -> Path:
    d = Path(user_data_dir(APP_NAME))
    d.mkdir(parents=True, exist_ok=True)
    return d
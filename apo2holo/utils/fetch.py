import requests
import subprocess
import shutil
import urllib.request
from pathlib import Path

RELEASE_BASE = "https://github.com/solarflip/hydrogenase-apo2holo-annotations/releases/tag/v0.0.1-structure-models"

REPO = "solarflip/hydrogenase-apo2holo-annotations"

FILES = {
    "as_experiment_02.pkl": f"{RELEASE_BASE}/as_experiment_02.pkl",
    "fes_pocket_experiment_02.pkl": f"{RELEASE_BASE}/fes_pocket_experiment_02.pkl",
    "fes_type_experiment_02.pkl": f"{RELEASE_BASE}/fes_type_experiment_02.pkl",
    "structures.zip": f"{RELEASE_BASE}/structures.zip"
}
# assets that need the more robust gh-cli download path instead of a plain HTTP GET
LARGE_ASSETS = {"structures.zip", "as_experiment_02.pkl", "fes_pocket_experiment_02.pkl", "fes_type_experiment_02.pkl"}

RELEASE_TAG = "v0.0.1-structure-models"

def fetch_file(name: str, dest_dir: Path, force: bool = False) -> Path:
    dest = dest_dir / name
    if dest.exists() and not force:
        return dest

    dest_dir.mkdir(parents=True, exist_ok=True)

    if name in LARGE_ASSETS:
        _fetch_via_gh_cli(name, dest_dir)
    else:
        _fetch_via_urllib(name, dest)

    return dest


def _fetch_via_urllib(name: str, dest: Path):
    url = FILES[name]
    print(f"Downloading {name} from {url} ...")
    urllib.request.urlretrieve(url, dest)


def _fetch_via_gh_cli(name: str, dest_dir: Path):
    if shutil.which("gh") is None:
        raise RuntimeError(
            f"'{name}' requires the GitHub CLI ('gh') to download. "
            "Install it from https://cli.github.com"
        )

    # confirm gh is authenticated (needed if the repo is private; harmless check if public)
    auth_check = subprocess.run(["gh", "auth", "status"], capture_output=True, text=True)
    if auth_check.returncode != 0:
        raise RuntimeError(
            "gh CLI is not authenticated. Run `gh auth login` first."
        )

    print(f"Downloading {name} via gh CLI (tag: {RELEASE_TAG}, repo: {REPO}) ...")
    result = subprocess.run(
        [
            "gh", "release", "download", RELEASE_TAG,
            "-R", REPO,
            "-p", name,
            "-D", str(dest_dir),
            "--clobber",
        ],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"gh release download failed for {name}:\n{result.stderr}"
        )

    dest = dest_dir / name
    if not dest.exists():
        raise RuntimeError(f"gh release download reported success but {dest} is missing")
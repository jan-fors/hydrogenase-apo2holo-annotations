"""
Given a chosen model each test data sample is predicted and results are saved
"""

import argparse
import os
from pathlib import Path
import subprocess
from src.main import main
from datetime import datetime

def predict(directory : Path, f_radius : float, model : Path, search_type : str, tmp : Path):
    """
    
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    model_name = model.stem
    runname = timestamp+ "__"+ model.parent.name +"__"+ model_name +"__"+search_type+"__"+str(f_radius)

    for subset in os.listdir(directory):
        print("annotate", subset, "...")
        subset_path = directory / Path(str(subset))
        db_dir = subset_path / Path('db')
        structure_db_path = db_dir / Path('structureDB') / Path('structureDB')
        fingerprint_db_dir_path = db_dir / Path('fingerprintDB')

        test_dir = subset_path / Path('test')

        
        model_absolut_path = fingerprint_db_dir_path / model

        pred_dir = subset_path / Path('out') / Path(runname)


        for test_structure in os.listdir(test_dir):
            try:
                print("Predicting test structure", test_structure)
                test_structure_path = test_dir / Path(test_structure)
                output = pred_dir / Path(test_structure.split(".")[0])
                os.makedirs(output, exist_ok=True)
                main(input_path=test_structure_path,
                    out=output,
                    tmp=tmp,
                    boltz=False,
                    plot=True, 
                    structure_db_path=structure_db_path,
                    fingerprint_db_path=fingerprint_db_dir_path,
                    search_type=search_type,
                    model=model_absolut_path,
                    f_radius=f_radius)
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    #parser.add_argument("name")
    parser.add_argument("--f_radius", type=float, default=5.0)
    parser.add_argument("--model", type=Path, help="Path inside the fingeprintDB directory to the model that should be used")
    parser.add_argument("--search_type", type=str, choices=["sum", "mc"], default="sum")
    parser.add_argument("--tmp", type=Path, default=Path("tmp"))

    args = parser.parse_args()

    predict(args.dir, args.f_radius, args.model, args.search_type, args.tmp)
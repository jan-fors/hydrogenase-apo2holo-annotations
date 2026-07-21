import os
from pathlib import Path
from src.utils.constants import STRUCTURE_DB, STRUCTURE_DIR, FOLDSEEK_OUT_FORMAT, TMP
from src.io.printl import printl
from src.utils.protein.get_chains import get_chains
from src.utils.protein.extract_chain import extract_chain
from src.filter.apply_blacklist import apply_blacklist_to_input_structures
from Bio import PDB
import subprocess
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed


class StructureDB:
    def __init__(
        self,
        input_data_dir_path: Path = None,
        output_data_dir_path: Path = None,
        structure_db_path: Path = None,
    ):
        """
        Initialize
        """
        self.input_data_dir_path = input_data_dir_path
        self.output_data_dir_path = output_data_dir_path
        self.structure_db_path = structure_db_path

    def load(
        self,
        input_data_dir_path: Path = None,
        output_data_dir_path: Path = None,
        structure_db_path: Path = None,
    ):
        """
        Load an existing database to be able to search it
        """
        self.input_data_dir_path = input_data_dir_path
        self.output_data_dir_path = output_data_dir_path
        self.structure_db_path = structure_db_path

    def search(self, query_path: Path, output_path: Path):
        """
        foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmpFolder
        """
        result_name = os.path.basename(query_path).split(".")[0] + "_ms_res"
        output_path = os.path.join(str(output_path), result_name)

        cmd = [
            "foldseek",
            "easy-search",
            str(query_path),
            str(self.structure_db_path),
            str(output_path),
            str(TMP),
            "--format-output",
            FOLDSEEK_OUT_FORMAT,
        ]

        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )

        return output_path

    def build_db(
        self,
        data_dir: Path,
        output_data_dir_path: Path = STRUCTURE_DIR,
        structure_db_path: Path = STRUCTURE_DB,
        threads: int = 1,
    ):
        """
        create a database
        """
        self._check_input_dir(data_dir)
        self._check_output_dir(output_data_dir_path)

        self.input_data_dir_path = data_dir
        self.output_data_dir_path = output_data_dir_path

        files = os.listdir(data_dir)

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(
                    self._process_file, file, data_dir, output_data_dir_path
                ): file
                for file in files
            }
            for future in as_completed(futures):
                file = futures[future]
                try:
                    future.result()
                except Exception as e:
                    printl(f"Failed on {file}: {e}")

        self._create_foldseek_db(output_data_dir_path, structure_db_path)

        self.structure_db_path = structure_db_path

    def _process_file(self, file: str, data_dir: Path, output_data_dir_path: Path):
        """ """
        file_path = os.path.join(data_dir, file)
        if not os.path.exists(file_path):
            printl(f"{file_path} does not exist.")

        # get chains
        chains = get_chains(file_path)

        # split each structure into subunits
        for chain in chains:
            result_structure = extract_chain(file_path, output_data_dir_path, chain)
            apply_blacklist_to_input_structures(result_structure)
            self._use_first_model(result_structure)

        # EXPERIMENT
        shutil.copy2(file_path, output_data_dir_path)

    def get_structure_path(self, name: str, chain_dir: Path):
        """ """
        structure_path = os.path.join(chain_dir, name + ".pdb")
        if not os.path.exists(structure_path):
            print(f"{structure_path} does not exist.")
            return None
        return Path(structure_path)

    def _check_input_dir(self, input_dir: Path):
        """ """
        # check if path exists
        if not os.path.exists(input_dir):
            raise ValueError(f"Path {input_dir} does not exist")

        # check if path contains pdb files
        # TODO

    def _check_output_dir(self, output_dir: str):
        """ """
        if not os.path.exists(output_dir):
            printl(f"{output_dir} does not exist. Creating it ...")
            os.makedirs(output_dir, exist_ok=True)

    def _use_first_model(self, structure_path: Path):
        """ """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", structure_path)

        # Nur das erste Modell auswählen (Index 0)
        first_model = structure[0]

        # Speichern
        io = PDB.PDBIO()
        io.set_structure(first_model)
        io.save(str(structure_path))

    def _create_foldseek_db(
        self, structure_dir_path: str, structure_db_path: str = STRUCTURE_DB
    ):
        """"""
        # create db
        cmd = ["foldseek", "createdb", str(structure_dir_path), structure_db_path]
        print(cmd)
        subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)

        # create index
        cmd = ["foldseek", "createindex", structure_db_path, "tmp"]
        subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)

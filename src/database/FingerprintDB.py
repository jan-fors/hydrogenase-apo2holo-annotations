import os

class FingerprintDB:
    def __init__(self, db_path: str):
        self.db_path = db_path

        if os.path.exists(db_path):
            print(f"Loading fingerprint database from {db_path} ...")
            self.db = self._load_db()
    
    def _load_db(self):
        """ """
        pass

    def search(self, fingerprint):
        """ """
        pass

    def build_db(self, data):
        """ """
        pass
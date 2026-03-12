from src.utils.constants import (
    FINGERPRINT_DB
)
import pandas as pd

def search_against_fingerprint_db(F : dict, fingerprint_db_path : str = FINGERPRINT_DB):
    """
    TODO load fingerprintDB once
    ALA	ARG	ASN	ASP	CYS	GLN	GLU	GLY	HIS	ILE	LEU	LYS	MET	PHE	PRO	SER	THR	TRP	TYR	VAL
    """
    DB = pd.read_csv(fingerprint_db_path, sep="\t")
    
   
    subset = DB[DB["ALA"] == F["ALA"]]
    subset = subset[subset["ARG"] == F["ARG"]]
    subset = subset[subset["ASN"] == F["ASN"]]
    subset = subset[subset["ASP"] == F["ASP"]]
    subset = subset[subset["CYS"] == F["CYS"]]
    subset = subset[subset["GLN"] == F["GLN"]]
    subset = subset[subset["GLU"] == F["GLU"]]
    subset = subset[subset["GLY"] == F["GLY"]]
    subset = subset[subset["HIS"] == F["HIS"]]
    subset = subset[subset["ILE"] == F["ILE"]]
    subset = subset[subset["LEU"] == F["LEU"]]
    subset = subset[subset["LYS"] == F["LYS"]]
    subset = subset[subset["MET"] == F["MET"]]
    subset = subset[subset["PHE"] == F["PHE"]]
    subset = subset[subset["PRO"] == F["PRO"]]
    subset = subset[subset["SER"] == F["SER"]]
    subset = subset[subset["THR"] == F["THR"]]
    subset = subset[subset["TRP"] == F["TRP"]]
    subset = subset[subset["TYR"] == F["TYR"]]
    subset = subset[subset["VAL"] == F["VAL"]]


    return subset["res_name"].to_list()
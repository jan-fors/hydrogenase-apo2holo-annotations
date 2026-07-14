"""
    constants.py
"""
# imports
from pathlib import Path
import numpy as np

# general
VERBOSE =           True

# paths
STRUCTURE_DB =      Path("db/structureDB/structureDB")
STRUCTURE_DIR =     Path("db/single_chains")

FINGERPRINT_DB =    Path("db/fingerprintDB") 
ACTIVE_SITE_MODEL = Path("db/fingerprintDB/as_mlp_classifier_model.pkl")
FES_MODEL =         Path("db/fingerprintDB/fes_mlp_classifier_model.pkl")

TMP =               Path("tmp")

# lists
COFACTOR_BLACKLIST = [
    "FE2", "MPD", "CL", "IMD", "CSD", "OCS", "PSW", "CA", "SE7", 
    "UOX","CO3", "SBY", "H2S", "HOH", "SO4", "MRD", "VK3", "DHI", 
    "TRS", "IMD", "CSO", "PEG", "LMT", "PO4", "NA", "LI", "MLA", "CSX", 
    "OXY", "MQ9", "GOL", "KR", "CMO", "MG", "FE1" # remove as well?
]

COFACTOR_WHITELIST = [
    "SF3", "F3S", "F4S", "SF4", "FSX"
]

ACTIVE_SITE_WHITELIST = [
    "NFV", "NWN", "NI", "3NI", "NFU", "NFO", "FNE", "FCO"
]

AA_ORDER = [
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL"
]

# search
NN_CLUSTERING_RADIUS =  4.0

FOLDSEEK_OUT_FORMAT =   "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"
FIDENT_THRESHOLD =      0.1
BITS_THRESHOLD =        500

# fingerprint
FINGERPRINT_RADIUS =    6.0

SOAP_N_MAX =            6
SOAP_L_MAX =            4
SOAP_SIGMA =            0.5
SOAP_RBF =              "gto"
SOAP_ALLOWED_SPECIES =  ["C", "N", "O", "S"]

AC_ELEMENTS =           ["C", "N", "O", "S"]
AC_SHELL_WIDTH =        1.0

# predict
SEARCH_TYPE =           "mlp_sum"   #"absolut", "logreg_sum", "svm_sum", "svm_mc", "logreg_mc", "mlp_mc", "mlp_sum"

# build
AMOUNT_ARTIFICAL_SAMPLES =      12
MIN_DIST_ARTIFICAL_SAMPLES =    6.0

AUGMENTATION_RADIUS =           4.0
N_AUGMENTATIONS =               12

# models TODO write one active site and one fes model later
MODEL = ""

# logreg params
LOGREG_MAX_ITER = 10000
LOGREG_SOLVER = 'saga'
LOGREG_CLASS_WEIGHT = None
LOGREG_C = np.float64(30.09558248946819)
LOGREG_L1_RATIO = np.float64(0.2665616636524313)
LOGREG_TOL = np.float64(2.52725657659196e-06)
LOGREG_FIT_INTERCEPT = True

# svm params
SVC_MAX_ITER=10000
SVC_KERNEL="rbf"
SVC_DEGREE=3
SVC_GAMMA="scale"
SVC_SHRINKING=True
SVC_TOLERANCE=0.002
SVC_C=1.0

# mlpclassifier params
MLP_ACTIVATION = 'relu'
MLP_SOLVER = "adam"
MLP_ALPHA = np.float64(3.695879010047527e-05)
MLP_MAX_ITER = 10000
MLP_BATCH_SIZE = "auto"
MLP_BETA_1 = np.float64(0.9181621909301618)
MLP_BETA_2 = np.float64(0.9835234744946176)
MLP_EARLY_STOPPING = False
MLP_HIDDEN_LAYERS = (500,250)
MLP_LEARNING_RATE = 'invscaling'
MLP_LEARNING_RATE_INIT = np.float64(0.000976339545835049)
MLP_MOMENTUM = np.float64(0.5492048991630394)
MLP_TOL = np.float64(0.00015392565509282595)
MLP_VALIDATION_FRACTION = np.float64(0.2512210788984352)

# randomforest params
RANF_BOOTSTRAP = True
RANF_CLASS_WEIGHT = None
RANF_MAX_DEPTH = 50
RANF_MAX_FEATURES = 'sqrt'
RANF_MAX_SAMPLES = np.float64(0.9569315944661306)
RANF_MIN_SAMPLES_LEAF = 1
RANF_MIN_SAMPLE_SPLIT = 5
RANF_N_ESTIMATORS = 110
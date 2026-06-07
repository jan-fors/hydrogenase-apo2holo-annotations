from pathlib import Path


VERBOSE=True

COFACTOR_BLACKLIST = [
    "FE2", "MPD", "CL", "IMD", "CSD", "OCS", "PSW", "CA", "SE7", 
    "UOX","CO3", "SBY", "H2S", "HOH", "SO4", "MRD", "VK3", "DHI", 
    "TRS", "IMD", "CSO", "PEG", "LMT", "PO4", "NA", "LI", "MLA", "CSX", 
    "OXY", "MQ9", "GOL", "KR", "CMO", "MG", "FE1" # remove as well?
]

COFACTOR_WHITELIST = [
    "SF3", "NFV", "F3S", "NWN", "NI", "3NI", "F4S", "SF4", "NFU", "NFO", "FNE", "FCO", "FSX"
]

AA_ORDER = [
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL"
]

STRUCTURE_DIR=Path("db/single_chains")#"dat/single_chains"#2026-03-06-raw_structures"

NN_CLUSTERING_RADIUS = 4.0

STRUCTURE_DB=Path("db/structureDB/structureDB")
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"

FIDENT_THRESHOLD=0.1
BITS_THRESHOLD=250

FINGERPRINT_RADIUS = 5
FINGERPRINT_DB=Path("db/fingerprintDB")
CLASS_LABEL_ORDER = []
MODEL = ""

TMP="tmp"

SEARCH_TYPE = "logreg_sum"   #"absolut", "logreg_sum", "svm_sum", "svm_mc", "logreg_mc", "mlp_mc", "mlp_sum"

# logreg params
LOGREG_MAX_ITER = 1000
LOGREG_SOLVER = 'lbfgs'
LOGREG_CLASS_WEIGHT = "balanced"
LOGREG_C = 0.055
LOGREG_L1_RATIO = 0.34
LOGREG_TOL = 6.93e-05
LOGREG_FIT_INTERCEPT = False

# svm params
SVC_MAX_ITER=1000
SVC_KERNEL="rbf"
SVC_DEGREE=3
SVC_GAMMA="scale"
SVC_SHRINKING=True
SVC_TOLERANCE=0.002
SVC_C=1.0

# mlpclassifier params
MLP_ACTIVATION = "relu"
MLP_SOLVER = "adam"
MLP_ALPHA = 1e-5
MLP_MAX_ITER = 5000
MLP_BATCH_SIZE = "auto"
MLP_BETA_1 = 0.91
MLP_BETA_2 = 0.99
MLP_EARLY_STOPPING = False
MLP_HIDDEN_LAYERS = (20,20,10,5)
MLP_LEARNING_RATE = "adaptive"
MLP_LEARNING_RATE_INIT = 0.004
MLP_MOMENTUM = 0.52
MLP_TOL = 0.00017
MLP_VALIDATION_FRACTION = 0.1158

# randomforest params
RANF_BOOTSTRAP = True
RANF_CLASS_WEIGHT = "balanced"
RANF_MAX_DEPTH = None
RANF_MAX_FEATURES = "log2"
RANF_MAX_SAMPLES = 0.91
RANF_MIN_SAMPLES_LEAF = 2
RANF_MIN_SAMPLE_SPLIT = 6
RANF_N_ESTIMATORS = 638
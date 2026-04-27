VERBOSE=True

COFACTOR_BLACKLIST = [
    "1C1O", "FE2", "MPD","CL", "IMD", "CSD", "OCS", "PSW", "CA", "SE7", 
    "UOX","CO3", "SBY", "CL", "H2S", "HOH", "SO4", "MRD", "VK3", "DHI", 
    "TRS", "IMD", "CSO", "PEG", "LMT", "PO4", "NA", "LI", "MLA", "CSX", 
    "OXY", "H2S", "MQ9", "GOL", "KR", "CMO", "MG" # remove as well?
]

COFACTOR_WHITELIST = [
    "CMO", "SF3", "NFV", "F3S", "NWN", "NI", "3NI", "F4S", "SF4", "NFR", "FE2", "NFU", "NFO", "FNE", "FCO", "FSX", "CYN"
]

AA_ORDER = [
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL"
]

STRUCTURE_DIR="dat/single_chains"#"dat/single_chains"#2026-03-06-raw_structures"

STRUCTURE_DB="db/structureDB/structureDB"
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"

FIDENT_THRESHOLD=0.1
BITS_THRESHOLD=500

FINGERPRINT_RADIUS = 5
FINGERPRINT_DB="db/fingerprintDB"
MODEL=""

TMP="tmp"

SEARCH_TYPE = "svm_sum"   #"absolut","logreg_sum", "svm_sum", "svm_mc", "logreg_mc", "nn_mc", "nn_sum"

# logreg params
MAX_ITER = 1000
SOLVER = 'lbfgs'
CLASS_WEIGHT = "balanced"

# svm params
SVC_MAX_ITER=1000
SVC_KERNEL="rbf"
SVC_DEGREE=3
SVC_GAMMA="scale"
SVC_SHRINKING=True
SVC_TOLERANCE=0.002
SVC_C=1.0

# mlpclassifier params
MLP_SOLVER = "adam"
MLP_ALPHA = 1e-5
MLP_MAX_ITER = 5000
MLP_HIDDEN_LAYERS = (20,20,10,5)

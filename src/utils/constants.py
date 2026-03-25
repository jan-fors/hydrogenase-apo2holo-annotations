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

STRUCTURE_DIR="dat/single_chains"#"dat/single_chains"#2026-03-06-raw_structures"

STRUCTURE_DB="db/structureDB/structureDB"
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"

FIDENT_THRESHOLD=0.1
BITS_THRESHOLD=500

FINGERPRINT_RADIUS = 5 # is too small
FINGERPRINT_DB="db/fingerprintDB"
MODEL=""

TMP="tmp"

MAX_ITER = 400
SOLVER = 'lbfgs'
SEARCH_TYPE = "absolut"#"absolut"#

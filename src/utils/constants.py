VERBOSE=True

COFACTOR_BLACKLIST = [
    "FE2", "MPD","CL", "IMD", "CSD", "OCS", "PSW", "CA", "SE7", "UOX","CO3", "SBY", "CL", "H2S", "HOH", "SO4", "MRD", "VK3", "DHI", "TRS", "IMD", "CSO", "PEG", "LMT", "PO4", "NA", "LI", "MLA", "CSX", "OXY", "H2S", "MQ9", "GOL", "KR", "MG" # remove as well?
]

STRUCTURE_DIR="dat/chains"

STRUCTURE_DB="db/structureDB"
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"

FIDENT_THRESHOLD=0.1
BITS_THRESHOLD=500

FINGERPRINT_RADIUS = 5
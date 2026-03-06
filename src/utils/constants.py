VERBOSE=True

COFACTOR_BLACKLIST = [
    "FE2", "MPD","CL", "IMD", "CSD", "OCS", "PSW", "CA", "SE7", "UOX","CO3", "SBY", "CL", "H2S", "HOH", "SO4", "MRD", "VK3", "DHI", "TRS", "IMD", "CSO", "PEG", "LMT", "PO4", "NA", "LI", "MLA", "CSX", "OXY", "H2S", "MQ9", "GOL", "KR", "MG" # remove as well?
]

COFACTOR_WHITELIST = [
    "CMO", "SF3", "NFV", "F3S", "NWN", "NI", "3NI", "F4S", "SF4", "NFR", "FE2", "NFU", "NFO", "FNE", "FCO", "FSX", "CYN"
]

STRUCTURE_DIR="dat/single_chains"

STRUCTURE_DB="db/structureDB"
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"

FIDENT_THRESHOLD=0.1
BITS_THRESHOLD=500

FINGERPRINT_RADIUS = 5
FINGERPRINT_DB="db/fingerprintDB/fingerprint_db.tsv"
VERBOSE=True

COFACTOR_BLACKLIST = [
    "HOH",
    "GOL",
    "H",
    "O",
    "MG" # remove as well?
]

STRUCTURE_DIR="dat/chains"

STRUCTURE_DB="db/structureDB"
FOLDSEEK_OUT_FORMAT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,u,t"
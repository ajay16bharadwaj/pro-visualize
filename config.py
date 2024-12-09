# config.py

# Default column names for analysis
DEFAULT_COL_DEPLABEL = 'Label2'       # Label column
DEFAULT_COL_DEPSIGNIF = 'Imputed.FDR'  # Significance column
DEFAULT_GENE_NAME_COLUMN = 'Gene Name' # Gene Name Column from UniProt Annotation

# Default thresholds
DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_LOG2FC_THRESHOLD = 0.6

# Color mapping for plots
COLOR_MAP = {
    'Down-regulated': '#0072B2',
    'Significant, no change': '#F0E442',
    'Up-regulated': '#D55E00',
    'Non-significant': '#999999'
}

# Supported organisms
ORGANISM_DICT = {
    "human": "hsapiens",
    "mouse": "mmusculus",
    "rat": "rnorvegicus",
    "zebrafish": "drerio",
    "fruit fly": "dmelanogaster",
    "worm": "celegans",
    "yeast": "scerevisiae",
    "arabidopsis": "athaliana",
    "pig": "sscrofa",
    "cow": "btaurus",
    "chicken": "ggallus"
}

# File upload settings
SUPPORTED_FILE_FORMATS = ['.csv', '.tsv']
DEFAULT_SEPARATOR = '\t'
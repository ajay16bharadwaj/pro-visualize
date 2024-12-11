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

# Centralized figure defaults
DEFAULT_VOLCANO_CONFIG = {
    "title": "Volcano Plot",
    "x_label": "Log2 Fold Change",
    "y_label": "-Log10 FDR",
    "label_size": 15,
    "colors": {
        "Up-regulated": "#D55E00",
        "Down-regulated": "#0072B2",
        "Non-significant": "#999999",
        "Significant, no change": "#F0E442",
    },
    "threshold_lines": True,
}

# heatmap_config.py

DEFAULT_HEATMAP_CONFIG = {
    "title": "Heatmap",
    "x_label": "Samples",
    "y_label": "Proteins",
    "color_gradient": "viridis",
    "height": 600,
    "width": 600,
    "title_font_size": 16,
    "x_label_font_size": 14,
    "y_label_font_size": 14,
}
# config.py
import seaborn as sns

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

DEFAULT_VIOLIN_CONFIG = {
    "title": "Violin Plot",
    "label_size": 15,
    "colors": {},
}

DEFAULT_PROTEIN_PER_SAMPLE_CONFIG = {
    "title": "Proteins Per Sample",
    "x_label": "Sample Name",
    "y_label": "Number of Proteins",
    "label_font_size": 15,
    "title_font_size": 20,
    "colors": {},  # Colors will be populated dynamically based on group
}

DEFAULT_GLOBAL_CONFIG = {
    "group_colors": {},
    "default_col_label": DEFAULT_COL_DEPLABEL,
    "default_col_significance": DEFAULT_COL_DEPSIGNIF,
    "default_gene_name_column": DEFAULT_GENE_NAME_COLUMN,
}

DEFAULT_DENSITY_CONFIG = {
    "palette": sns.color_palette("husl", 10),  # Default color palette
    "height": 4,                               # Height of each plot
    "col_wrap": 3,                             # Number of columns in FacetGrid
    "alpha": 0.8,                              # Transparency of density lines
    "linewidth": 1.2                           # Line thickness of density plots
}

DEFAULT_PCA_BY_ANNOTATION_CONFIG = {
    "title": "PCA Plot - Colored by Group",
    "x_label": "PC1",
    "y_label": "PC2",
    "marker_size": 8,
    "marker_symbol": "circle",
    "jitter": 0.3,  # Add jitter to avoid overlap
    "width": 800,
    "height": 600,
}

DEFAULT_SANKEY_CONFIG = {
    "colors": {
        "Upregulated": "#00FF00",  # Green in hex
        "Downregulated": "#FF0000",  # Red in hex
        "Non-significant": "#D3D3D3"  # Light gray in hex
    },
    "pad": 20,
    "thickness": 20,
    "scale_non_significant": 0.3,
    "filter_threshold": 1,
    "font_size": 12,
    "title_font_size": 14,
    "title": "Sankey Diagram of Protein Transitions"
}

DEFAULT_UMAP_CONFIG = {
    "title": "UMAP Projection",
    "marker_size": 10,
    "width": 800,
    "height": 600,
    "colors": {},  # Dynamically populated from global config
}

DEFAULT_TSNE_CONFIG = {
    "title": "t-SNE Projection",
    "marker_size": 8,
    "width": 800,
    "height": 600,
    "colors": {},  # Dynamically populated from global config
    "show_labels": True,  # Toggle option for sample labels
}
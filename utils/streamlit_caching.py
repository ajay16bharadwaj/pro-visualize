import streamlit as st
import pandas as pd
from config import (
    DEFAULT_COL_DEPLABEL,
    DEFAULT_COL_DEPSIGNIF,
    COLOR_MAP,
    ORGANISM_DICT,
    DEFAULT_FDR_THRESHOLD,
    DEFAULT_LOG2FC_THRESHOLD,
    DEFAULT_SEPARATOR
)

######################################
# functions - used for caching 
######################################
@st.cache_data
def load_data(path):
    df = pd.read_csv(path,sep = DEFAULT_SEPARATOR)
    return df

@st.cache_data 
def cached_get_all_enrichment(_vis, protein_list, organism):
    """
    Cached version of the enrichment analysis method.
    This function calls the method from the `vis` object and caches the result.
    
    Args:
        vis: Instance of the class containing the `get_all_enrichment` method.
        protein_list (list): List of proteins for enrichment analysis.
        organism (str): The organism for the analysis (e.g., 'human').

    Returns:
        pd.DataFrame: The complete enrichment DataFrame.
        dict: A dictionary of DataFrames for different sources.
    """
    enrichment_df, source_dict =  _vis.get_all_enrichment(protein_list, organism)
    return enrichment_df, source_dict
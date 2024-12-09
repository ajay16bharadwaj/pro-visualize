from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd

#Volcano Plot - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("Volcano Plot")
@validate_inputs(analysis_status=True, protein_status=True, annotation_status=True)
def render_volcano_plot(vis, figures_dict, analysis_status, protein_status, annotation_status, **kwargs):
    """Render Volcano Plot Tab and update figures_dict."""
    condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL)
    volcano_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, vis.dep_info['log2FC'].max(), value=0.6, key='volcano_fc') # type: ignore
    volcano_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='volcano_fdr')
    volcano_comparison_input = st.selectbox('Which comparison', condition_groups, index=0, key="volcano_comparison_input" )
    
    # Generate figure
    fig = vis.plot_volcano(condition_groups[volcano_comparison_input], volcano_log2fc_threshold, volcano_fdr_threshold)
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    figures_dict["volcano_plot"] = fig


from matplotlib import pyplot as plt
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from io import BytesIO
import numpy as np
import dash_bio
from sklearn.preprocessing import MinMaxScaler
from compute_functions import create_download_button
from tabs.clustering import render_clustering
from tabs.functional_annotation import render_functional_annotation
from tabs.quantifications import render_quantification_plots
from tabs.uniprot_annotation import render_get_uniprot_annotation
from tabs.venn_tab import render_venn_tab
from tabs.violin_plot import render_violin_plot
from tabs.volcano_plot import render_volcano_plot
from tabs.heatmap import render_heatmap
from visualization import ProteinVisualization
from config import (
    DEFAULT_COL_DEPLABEL,
    DEFAULT_COL_DEPSIGNIF,
    DEFAULT_GENE_NAME_COLUMN,
    COLOR_MAP,
    ORGANISM_DICT,
    DEFAULT_FDR_THRESHOLD,
    DEFAULT_LOG2FC_THRESHOLD,
    DEFAULT_SEPARATOR,
    DEFAULT_GLOBAL_CONFIG,
)



#######################################
# PAGE SETUP
#######################################

st.set_page_config(page_title="ProEpic Report", page_icon=":bar_chart:", layout="wide")

st.title("Pro-Visualize")
st.markdown("_Prototype v0.4.1_")

analysis_status = False
protein_level_status = False
annotation_status = False

with st.sidebar:

    st.header("Inputs")
    analysis_upload = st.file_uploader("Analysis Input")

    if analysis_upload is None:
        st.info(" Upload Differentially expressed file through config", icon="ℹ️")

    protein_level_upload = st.file_uploader("Protein Level Input")

    if protein_level_upload is None:
        st.info( " Upload the Protein Level input through config", icon="ℹ️")

    annotation_file_upload = st.file_uploader("Upload your metadata/Annotation file")

    if annotation_file_upload is None:
        st.info( " Upload the metadata/Annotation file through config", icon="ℹ️")
        


#######################################
# DATA LOADING
#######################################

# Initialize dictionary to store the figures
figures_dict = {}
#intialize class 
vis = ProteinVisualization() 

with st.container():

    preview_col1, preview_col2, preview_col3 = st.columns(3)

    with preview_col1:
        with st.expander("Analysis Input Preview"):
            if analysis_upload is not None:
                #analysis_df = load_data(analysis_upload) # type: ignore
                vis.load_dep_data(analysis_upload) # type: ignore
                #vis.dep_info = analysis_df.copy()
                #vis.analysis_with_annotation = analysis_df.copy()
                st.dataframe(vis.dep_info)
                analysis_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")
    with preview_col2:
        with st.expander("Annotation Preview"):
            if annotation_file_upload is not None:
                #annotation_file_upload = load_data(annotation_file_upload)
                vis.load_annotation(annotation_file_upload, "Level3", "attribute_ExperimentalGroup" ) # type: ignore #parameter values later on. 
                st.dataframe(vis.annotation_info)
                annotation_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")

    with preview_col3:
        with st.expander("Protein Input Preview"):
            if protein_level_upload is not None:
                #protein_level_upload = load_data(protein_level_upload)
                vis.load_protein_data(protein_level_upload) # type: ignore
                #upload pt level with uniport annotation - to save time. 
                vis.protein_data_annotated = vis.protein_data.copy()
                st.dataframe(vis.protein_data_annotated)
                protein_level_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")
 
# Default settings
global_config = DEFAULT_GLOBAL_CONFIG.copy()
if "global_config" not in st.session_state:
            st.session_state.global_config = DEFAULT_GLOBAL_CONFIG.copy()

with st.expander("⚙️ Global Configuration", expanded=False):
    st.header("Global Configuration")
    

    # Default settings
    global_config["default_col_label"] = st.text_input(
        "Default Label Column", 
        global_config["default_col_label"], 
        key="global_col_label"
    )
    global_config["default_col_significance"] = st.text_input(
        "Default Significance Column", 
        global_config["default_col_significance"], 
        key="global_col_significance"
    )
    global_config["default_gene_name_column"] = st.text_input(
        "Default Gene Name Column", 
        global_config["default_gene_name_column"], 
        key="global_gene_name_column"
    )

    # Customize group colors
    st.subheader("Group Colors")
    unique_groups = vis.annotation_info['Group'].unique() if vis.annotation_info is not None else []
    default_colors = px.colors.qualitative.Plotly[:len(unique_groups)]
    for i, group in enumerate(unique_groups):
        default_color = default_colors[i] if i < len(default_colors) else f"#%06x" % (0xFFFFFF & hash(group))
        global_config["group_colors"][group] = st.color_picker(
            f"Color for {group}",
            default_color,
            key=f"global_group_color_{group}"
        )

    # Apply and Reset Buttons
    apply_col, reset_col = st.columns([1, 1])
    with apply_col:
        if st.button("Apply Changes", key="global_apply_changes"):
            st.session_state.global_config = global_config.copy()
            st.toast("Global configuration updated!", icon="✅")
    with reset_col:
        if st.button("Reset to Default", key="global_reset_default"):
            st.session_state.global_config = DEFAULT_GLOBAL_CONFIG.copy()
            st.session_state.global_config['colors'] = default_colors



#initializing tabs
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(["Volcano Plot", "Heat Maps", "Violin Plot", "Quantification", "Clustering", "Venn Diagram", "Functional Analysis and Biological Annotations", "Get Uniprot annotations"])

#volcano plot - to set interactive FDR and FC thresholds
with tab1:
    render_volcano_plot(vis, figures_dict, analysis_status=analysis_status, protein_status=protein_level_status, annotation_status=annotation_status)
      

#heatmaps
with tab2:
    render_heatmap(vis, figures_dict, analysis_status=analysis_status, protein_status=protein_level_status, annotation_status=annotation_status)

#violin plots
with tab3:
    #st.write("Violin Plots")
    render_violin_plot(vis, figures_dict, analysis_status=analysis_status, protein_status=protein_level_status, annotation_status=annotation_status, global_config=st.session_state.global_config)

#quantification plots 
with tab4:
    #st.write("quantification")
    render_quantification_plots(vis, figures_dict, protein_status=protein_level_status, annotation_status=annotation_status, global_config=st.session_state.global_config)
    
#clustering
with tab5:
    render_clustering(vis, figures_dict, protein_status=protein_level_status, annotation_status=annotation_status, global_config=st.session_state.global_config)

#venn diagram
with tab6:
    #st.write("Venn Diagram")
    render_venn_tab(vis, figures_dict, protein_status=protein_level_status, annotation_status=annotation_status)

#functional annotation
with tab7:
    #st.write("Functional analysis and biological annotation")
    render_functional_annotation(vis, figures_dict, analysis_status=analysis_status, protein_status=protein_level_status, annotation_status=annotation_status)

#uniprot annotation
with tab8:
    st.write("Uniprot Annotations")
    render_get_uniprot_annotation(vis)

with st.sidebar:
    create_download_button(figures_dict)

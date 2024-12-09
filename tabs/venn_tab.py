from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd
from io import BytesIO
from utils.helpers import dataframe_with_selections

#quantification plots - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("venn")
@validate_inputs(protein_status=True, annotation_status=True)
def render_venn_tab(vis, figures_dict, protein_status, annotation_status, **kwargs):
    """Render venn diagrams and update figures_dict."""
    #select box for choosing custom groups
    #venn_custom_group_select_checkbox = st.checkbox(' Choose custom groups ', key='venn_group_select')
    venn_tab1, venn_tab2 = st.tabs(['Venn Diagram for Proteins Identified', 'Venn Diagram across comparisons'])
    with venn_tab1: 
        #st.write("For proteins identified")
        
        group_map = vis.annotation_info.set_index("SampleName")["Group"].to_dict()
        grouped_samples = vis.annotation_info.groupby("Group")["SampleName"].apply(list)
        all_groups = list(grouped_samples.keys())
        selected_groups = st.multiselect("Select groups to include", all_groups, default=all_groups, key='venn_tab1_select')

        # Filter to include only the selected groups.
        selected_samples = {group: samples for group, samples in grouped_samples.items() if group in selected_groups}

        # Step 4: Create a dictionary with protein sets for each selected group.
        protein_grouped_list = {}
        for group, samples in selected_samples.items():
            # Get a subset of protein_data that includes only the group's samples.
            group_proteins = vis.protein_data[["ProteinIds"] + samples]
            
            # Filter out proteins with NaN values across all samples for that group.
            identified_proteins = group_proteins.dropna(subset=samples, how="all")["ProteinIds"]
            
            # Store the set of identified proteins for this group.
            protein_grouped_list[group] = set(identified_proteins)

        venn_diagram = vis.plot_venn(protein_grouped_list)
        # Save the figure to a BytesIO object.
        img_bytes = BytesIO()
        venn_diagram.savefig(img_bytes, format='png', bbox_inches='tight')
        img_bytes.seek(0)

        # Use st.image to display the image with a specified width.
        st.image(img_bytes, caption='Venn Diagram', use_column_width=False, width=600)
    
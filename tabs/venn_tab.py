from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_SANKEY_CONFIG 
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
    venn_tab1, venn_tab2, venn_tab3 = st.tabs(['Venn Diagram for Proteins Identified', 'Sankey Plot', 'upSet Plot'])
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

    with venn_tab2: 
        st.subheader("Sankey Diagram")

        # Threshold sliders
        log2fc_threshold = st.slider(
            "Log2 Fold Change Threshold", 
            0.0, 5.0, 0.6, 0.1, key='sankey_logfc_threshold'
        )
        fdr_threshold = st.slider(
            "FDR Threshold", 
            0.0, 1.0, 0.05, 0.01, key='sankey_fdr_threshold'
        )
        

        # Initialize session state for Sankey config
        if "sankey_config" not in st.session_state:
            st.session_state.sankey_config = DEFAULT_SANKEY_CONFIG.copy()

        # Customization options
        with st.expander("Customize Sankey Diagram", expanded=False):
            pad = st.slider(
                "Node Padding", 
                0, 50, st.session_state.sankey_config.get("pad", 20), 1, key='sankey_node_padding'
            )
            thickness = st.slider(
                "Node Thickness", 
                10, 50, st.session_state.sankey_config.get("thickness", 20), 1, key='sankey_node_thickness'
            )
            scale_non_significant = st.slider(
            "Scale Non-significant Flows", 
            0.0, 1.0, st.session_state.sankey_config.get("scale_non_significant", 0.3), 0.1
            )
            filter_threshold = st.number_input(
                "Minimum Protein Count for Links", 
                min_value=1, value=st.session_state.sankey_config.get("filter_threshold", 1), step=1
            )
            font_size = st.slider(
                "Font Size", 
                8, 20, st.session_state.sankey_config.get("font_size", 12), 1, key='sankey_font_size'
            )
            title_font_size = st.slider(
                "Title Font Size", 
                10, 30, st.session_state.sankey_config.get("title_font_size", 14), 1, key='sankey_title_font_size'
            )

            # Add color pickers for categories
            upregulated_color = st.color_picker(
                "Upregulated Color", 
                st.session_state.sankey_config["colors"]["Upregulated"], key='sankey_upregulated_color'
            )
            downregulated_color = st.color_picker(
                "Downregulated Color", 
                st.session_state.sankey_config["colors"]["Downregulated"], key='sankey_downregulated_color'
            )
            non_significant_color = st.color_picker(
                "Non-significant Color", 
                st.session_state.sankey_config["colors"]["Non-significant"], key='sankey_non_significant_color'
            )
        
        # Buttons for applying and resetting
        apply_col, reset_col = st.columns(2)
        
        with apply_col:
            if st.button("Apply Changes", key="sankey_apply_changes"):
                st.session_state.sankey_config.update({
                    "pad": pad,
                    "thickness": thickness,
                    "font_size": font_size,
                    "title_font_size": title_font_size,
                    "colors": {
                        "Upregulated": upregulated_color,
                        "Downregulated": downregulated_color,
                        "Non-significant": non_significant_color
                    }
                })
                st.toast("Sankey diagram customization applied!", icon="✅")

        with reset_col:
            if st.button("Reset to Defaults", key="sankey_reset_defaults"):
                st.session_state.sankey_config = DEFAULT_SANKEY_CONFIG.copy()
                st.toast("Sankey diagram reset to default settings!", icon="🔄")
        
        # Generate Sankey plot with updated config
        sankey_plot = vis.plot_sankey(
            log2fc_threshold=log2fc_threshold,
            fdr_threshold=fdr_threshold,
            scale_non_significant=scale_non_significant,
            filter_threshold=filter_threshold,
            config=st.session_state.sankey_config
        )
        st.plotly_chart(sankey_plot, use_container_width=True)

    with venn_tab3:
        
        available_groups = vis.dep_info[DEFAULT_COL_DEPLABEL].unique().tolist()

        # Add a multiselect box for users to select groups
        selected_groups = st.multiselect(
            "Select groups to include in the UpSet plot",
            available_groups,
            default=available_groups,  # Pre-select all groups by default
            key="upset_group_selector"
        )

        # Generate and display the UpSet Plot
        with st.spinner("Generating UpSet plot..."):
            upset_image = vis.plot_upset(selected_groups=selected_groups)
            st.image(upset_image, caption="UpSet Plot of Protein Intersections", use_column_width=True)
    
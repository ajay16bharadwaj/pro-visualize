from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_PROTEIN_PER_SAMPLE_CONFIG, DEFAULT_PROTEIN_INTENSITY_DENSITY_CONFIG
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

#quantification plots - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("Quantification")
@validate_inputs(protein_status=True, annotation_status=True)
def render_quantification_plots(vis, figures_dict, protein_status, annotation_status, global_config, **kwargs):
    """Render Quantification plots Tab and update figures_dict."""
    quant_tab1, quant_tab2, quant_tab3, quant_tab4, quant_tab5 = st.tabs(['Protein Per Sample', 'Protein Overlap', 'Protein Intensity Density', 'Correlation Matrix', 'Protein Rank Order'])
        
    with quant_tab1:
        # Initialize session state for plot configuration
        if "protein_per_sample_config" not in st.session_state:
            st.session_state.protein_per_sample_config = DEFAULT_PROTEIN_PER_SAMPLE_CONFIG.copy()

        # Allow user to modify plot settings
        with st.expander("Modify Plot"):
            st.session_state.protein_per_sample_config["title"] = st.text_input(
                "Plot Title", 
                st.session_state.protein_per_sample_config["title"], 
                key="protein_per_sample_title"
            )
            st.session_state.protein_per_sample_config["x_label"] = st.text_input(
                "X-axis Label", 
                st.session_state.protein_per_sample_config["x_label"], 
                key="protein_per_sample_x_label"
            )
            st.session_state.protein_per_sample_config["y_label"] = st.text_input(
                "Y-axis Label", 
                st.session_state.protein_per_sample_config["y_label"], 
                key="protein_per_sample_y_label"
            )
            st.session_state.protein_per_sample_config["label_font_size"] = st.slider(
                "Label Font Size", 
                10, 30, 
                st.session_state.protein_per_sample_config["label_font_size"], 
                key="protein_per_sample_label_font_size"
            )
            st.session_state.protein_per_sample_config["title_font_size"] = st.slider(
                "Title Font Size", 
                10, 40, 
                st.session_state.protein_per_sample_config["title_font_size"], 
                key="protein_per_sample_title_font_size"
            )

            # Apply and Reset Buttons
            apply_col, reset_col = st.columns([1, 1])
            with apply_col:
                if st.button("Apply Changes", key="quant_apply_changes"):
                    st.toast("Quantification plot configuration updated!", icon="✅")
            with reset_col:
                if st.button("Reset to Defaults", key="quant_reset_defaults"):
                    st.session_state.protein_per_sample_config = DEFAULT_PROTEIN_PER_SAMPLE_CONFIG.copy()
                    st.toast("Quantification plot configuration reset to defaults!", icon="🔄")

        # Generate the plot using user-modified configuration and global_config for colors
        plot_proteins_per_sample = vis.plot_proteins_per_sample(
            group_column='Group',
            color_discrete_map=global_config["group_colors"],  # Colors from global_config
            config=st.session_state.protein_per_sample_config  # User-modified plot config
        )
        st.plotly_chart(plot_proteins_per_sample)
        figures_dict["proteins_per_sample"] = plot_proteins_per_sample

    with quant_tab2: 
        plot_protein_overlap = vis.plot_protein_overlap()
        st.plotly_chart(plot_protein_overlap)
        figures_dict["protein_overlap"] = plot_protein_overlap

    with quant_tab3:
        if "intensity_density_config" not in st.session_state:
            st.session_state.intensity_density_config = DEFAULT_PROTEIN_INTENSITY_DENSITY_CONFIG.copy()

        # Create a temporary structure to hold user inputs
        user_inputs = st.session_state.intensity_density_config.copy()

        # User inputs for customizing the plot
        with st.expander("Modify Density Plot"):
            user_inputs["title"] = st.text_input(
                "Plot Title", 
                user_inputs["title"], 
                key="density_plot_title"
            )
            user_inputs["x_label"] = st.text_input(
                "X-axis Label", 
                user_inputs["x_label"], 
                key="density_x_label"
            )
            user_inputs["y_label"] = st.text_input(
                "Y-axis Label", 
                user_inputs["y_label"], 
                key="density_y_label"
            )
            user_inputs["nbins"] = st.slider(
                "Number of Bins", 
                10, 200, 
                user_inputs["nbins"], 
                key="density_nbins"
            )
            user_inputs["height"] = st.slider(
                "Plot Height", 
                400, 1200, 
                user_inputs["height"], 
                key="density_height"
            )
            user_inputs["width"] = st.slider(
                "Plot Width", 
                400, 1200, 
                user_inputs["width"], 
                key="density_width"
            )

        # Apply and Reset buttons
        apply_col, reset_col = st.columns([1, 1])
        with apply_col:
            if st.button("Apply Changes", key="density_apply_changes"):
                # Save the user inputs into session state
                st.session_state.intensity_density_config = user_inputs.copy()
                st.toast("Density plot configuration updated!", icon="✅")
        with reset_col:
            if st.button("Reset to Default", key="density_reset_default"):
                # Reset to the default configuration
                st.session_state.intensity_density_config = DEFAULT_PROTEIN_INTENSITY_DENSITY_CONFIG.copy()
                st.toast("Density plot configuration reset to defaults!", icon="🔄")

        # Generate the plot using the updated configuration
        intensity_density_fig = vis.plot_intensity_density(
            config=st.session_state.intensity_density_config,
            color_discrete_map=global_config["group_colors"]
        )
        st.plotly_chart(intensity_density_fig)
        figures_dict["intensity_density"] = intensity_density_fig
        # plot_intensity_density = vis.plot_intensity_density()
        # st.plotly_chart(plot_intensity_density)
        # figures_dict["intensity_density"] = plot_intensity_density

    with quant_tab4:
        plot_correlation_matrix = vis.plot_correlation_matrix()
        st.plotly_chart(plot_correlation_matrix)
        figures_dict["correlation_matrix"] = plot_correlation_matrix
        

    with quant_tab5: 
        protein_highlight_select = st.checkbox('Choose Proteins to highlight', key='protein_rank_order_custom_select')
        #if certain proteins need to be highlighted 
        if protein_highlight_select:
            selection = dataframe_with_selections(vis.protein_data, "protein_rank_order_custom_df_select")
            with st.expander("Your selection"):
                st.write(selection)

            
            selected_proteins = list(selection['Protein'])
            protein_rank_order_plot = vis.plot_protein_rank_order(selected_proteins)
            st.plotly_chart(protein_rank_order_plot, use_container_width=True)
            figures_dict["protein_rank_order"] = protein_rank_order_plot

        else:
            protein_rank_order_plot = vis.plot_protein_rank_order()
            st.plotly_chart(protein_rank_order_plot, use_container_width=True)
            figures_dict["protein_rank_order"] = protein_rank_order_plot
    
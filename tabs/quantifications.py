from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

#quantification plots - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("Quantification")
@validate_inputs(protein_status=True, annotation_status=True)
def render_quantification_plots(vis, figures_dict, protein_status, annotation_status, **kwargs):
    """Render Quantification plots Tab and update figures_dict."""
    quant_tab1, quant_tab2, quant_tab3, quant_tab4, quant_tab5 = st.tabs(['Protein Per Sample', 'Protein Overlap', 'Protein Intensity Density', 'Correlation Matrix', 'Protein Rank Order'])
        
    with quant_tab1: 
        plot_proteins_per_sample =  vis.plot_proteins_per_sample()
        st.plotly_chart(plot_proteins_per_sample)
        figures_dict["proteins_per_sample"] = plot_proteins_per_sample

    with quant_tab2: 
        plot_protein_overlap = vis.plot_protein_overlap()
        st.plotly_chart(plot_protein_overlap)
        figures_dict["protein_overlap"] = plot_protein_overlap

    with quant_tab3:
        plot_intensity_density = vis.plot_intensity_density()
        st.plotly_chart(plot_intensity_density)
        figures_dict["intensity_density"] = plot_intensity_density

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
    
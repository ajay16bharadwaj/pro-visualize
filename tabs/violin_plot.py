from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

#violin plot - Tab Code. Need documentation for what the plot is? 
#default to top 10 DEP, and options FDR and FC values. Give the option to select custom proteins as well. 
@safe_tab_execution("Violin Plot")
@validate_inputs(analysis_status=True, protein_status=True, annotation_status=True)
def render_violin_plot(vis, figures_dict, analysis_status, protein_status, annotation_status, **kwargs):
    """Render Violin plots Tab and update figures_dict."""
    violin_info = vis.dep_info.copy() #type:ignore
    violin_condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL) # type: ignore
    violin_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, violin_info['log2FC'].max(), value=0.6, key='violin_fc') # type: ignore
    violin_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='violin_fdr')
    violin_comparison_input = st.selectbox('Which comparison', violin_condition_groups, index=0, key="violin_comparison_input" )

    #getting the top 10 differentially expressed proteins and plotting that 
    df = violin_condition_groups[violin_comparison_input]
    dep_list_df = df[(df['Imputed.FDR'] < violin_fdr_threshold) & ((df['log2FC'] < -(violin_log2fc_threshold)) | (df['log2FC'] > violin_log2fc_threshold)) ]

    #protein level preprocessing for violin plot
    pt_level_for_violin = vis.preprocess_for_violin_plot()

    with st.expander("Differentially expressed proteins"):
        st.dataframe(dep_list_df)

    custom_row_select = st.checkbox(' Choose custom entries ', key='violin_custom_select')

    if custom_row_select:
        subset_df = df[['Protein', 'log2FC', 'Imputed.FDR', 'Gene Name', 'Protein Description']]
        selection = dataframe_with_selections(subset_df, "violin_custom_df_select")
        with st.expander("Your selection"):
            st.write(selection)

        
        selected_df = pt_level_for_violin[pt_level_for_violin['Protein'].isin(list(selection['Protein']))]
        if (len(selected_df) <= 1):
            st.info("Please select at-least two inputs to continue")
        else:
            fig = vis.plot_violin_with_subplots(selected_df['Protein'].to_list())
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["violin_plot"] = fig

    else:
        dep_list_df['abs_log2FC'] = dep_list_df['log2FC'].abs()  # Create a new column with the absolute values of log2FC
        top10_by_fc = dep_list_df.sort_values(by='abs_log2FC', ascending=False).head(10)
        top_10_df = pt_level_for_violin[pt_level_for_violin['Protein'].isin(list(top10_by_fc['Protein'])[:10])]

        with st.expander("Top 10 DEP"):
            st.dataframe(top_10_df)

        if (len(top_10_df) == 0):
            st.info("Your selection has resulted in empty dataframe, please select different thresholds and try again")
        else:
            fig = vis.plot_violin_with_subplots(top_10_df['Protein'].to_list())
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["violin_plot"] = fig
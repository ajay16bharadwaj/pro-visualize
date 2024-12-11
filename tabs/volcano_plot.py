from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_VOLCANO_CONFIG 
import streamlit as st
import pandas as pd
from utils.helpers import get_user_volcano_config, dataframe_with_selections

#Volcano Plot - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("Volcano Plot")
@validate_inputs(analysis_status=True, protein_status=True, annotation_status=True)
def render_volcano_plot(vis, figures_dict, analysis_status, protein_status, annotation_status, **kwargs):
    """Render Volcano Plot Tab and update figures_dict."""
    condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL)
    # Initialize session state for user_config
    if "user_config" not in st.session_state:
        st.session_state.user_config = DEFAULT_VOLCANO_CONFIG.copy()


    selected_genes_for_volcano = []
    # Threshold and comparison inputs
    col1, col2 = st.columns(2)
    with col1:
        volcano_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, vis.dep_info['log2FC'].max(), 0.6)
    with col2:
        volcano_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05)
    volcano_comparison_input = st.selectbox('Which Comparison', condition_groups, index=0)

    # Collect user inputs directly
    with st.expander("Customize Volcano Plot"):
        title = st.text_input("Plot Title", st.session_state.user_config["title"])
        x_label = st.text_input("X-axis Label", st.session_state.user_config["x_label"])
        y_label = st.text_input("Y-axis Label", st.session_state.user_config["y_label"])
        label_size = st.slider(
            "Label Font Size", 
            10, 30, 
            st.session_state.user_config["label_size"]
        )
        
        # Organize color pickers into two columns
        color_col1, color_col2 = st.columns(2)
        with color_col1:
            upregulated_color = st.color_picker(
                "Up-regulated Color", 
                st.session_state.user_config["colors"]["Up-regulated"]
            )
            downregulated_color = st.color_picker(
                "Down-regulated Color", 
                st.session_state.user_config["colors"]["Down-regulated"]
            )
        with color_col2:
            nonsignificant_color = st.color_picker(
                "Non-significant Color", 
                st.session_state.user_config["colors"]["Non-significant"]
            )
            significant_no_change_color = st.color_picker(
                "Significant, no change Color", 
                st.session_state.user_config["colors"]["Significant, no change"]
            )
        
        threshold_lines = st.checkbox(
            "Show Threshold Lines", 
            st.session_state.user_config["threshold_lines"]
        )

    protein_highlight_select = st.checkbox('Choose Gene Names to highlight', key='volcano_plot_custom_select')
    #if certain proteins need to be highlighted 
    if protein_highlight_select:
        with st.expander('Choose Gene Names dropdown'):
            filtered_df = condition_groups[volcano_comparison_input]
            #getting the list of differentially expressed proteins for that comparison
            dep_list_df = filtered_df[(filtered_df['Imputed.FDR'] < volcano_fdr_threshold) & ((filtered_df['log2FC'] < -(volcano_log2fc_threshold)) | (filtered_df['log2FC'] > volcano_log2fc_threshold)) ]
            if dep_list_df.empty:
                st.warning("No differentially expressed proteins found for the selected comparison. Please select a different comparison.")
                st.stop()  # Stops the execution here, so the user has to select a valid comparison.
            
            selection = dataframe_with_selections(dep_list_df, "volcano_plot_custom_df_select")
            
            selected_genes_for_volcano = list(selection['Gene Name'])

    # Buttons for applying changes and resetting to defaults
    apply_col, reset_col = st.columns([1, 1])
    with apply_col:
        if st.button("Apply Changes", key='volcano_apply_changes'):
            st.session_state.user_config.update({
                "title": title,
                "x_label": x_label,
                "y_label": y_label,
                "label_size": label_size,
                "colors": {
                    "Up-regulated": upregulated_color,
                    "Down-regulated": downregulated_color,
                    "Non-significant": nonsignificant_color,
                    "Significant, no change": significant_no_change_color,
                },
                "threshold_lines": threshold_lines,
            })
            st.toast("Plot customization updated!", icon="✅")
            
    with reset_col:
        if st.button("Reset to Defaults", key='volcano_reset_defaults'):
            st.session_state.user_config = DEFAULT_VOLCANO_CONFIG.copy()
            st.toast("Reset to default settings!", icon="🔄")    

    # Generate the plot
    fig = vis.plot_volcano(condition_groups[volcano_comparison_input], volcano_log2fc_threshold, volcano_fdr_threshold, selected_genes_for_volcano, st.session_state.user_config)
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)


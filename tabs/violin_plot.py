from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_VIOLIN_CONFIG
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections
import plotly.express as px

@safe_tab_execution("Violin Plot")
@validate_inputs(analysis_status=True, protein_status=True, annotation_status=True)
def render_violin_plot(vis, figures_dict, analysis_status, protein_status, annotation_status, **kwargs):
    """Render Violin Plot Tab and update figures_dict."""
    # Preprocess data for plotting
    violin_info = vis.dep_info.copy()  # type: ignore
    violin_condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL)  # type: ignore
    
    # Initialize session state for user configurations
    if "violin_config" not in st.session_state:
        st.session_state.violin_config = DEFAULT_VIOLIN_CONFIG.copy()


    # Threshold and comparison inputs
    col1, col2 = st.columns(2)
    with col1:
        violin_log2fc_threshold = st.slider(
            'Log2 Fold Change Threshold', 
            0.0, violin_info['log2FC'].max(), 
            value=0.6, 
            key='violin_fc'
        )
    with col2:
        violin_fdr_threshold = st.slider(
            'Imputed FDR Threshold', 
            0.0, 1.0, 
            value=0.05, 
            key='violin_fdr'
        )
    violin_comparison_input = st.selectbox(
        'Which comparison', 
        violin_condition_groups, 
        index=0, 
        key="violin_comparison_input"
    )

    # Differentially expressed proteins
    df = violin_condition_groups[violin_comparison_input]
    dep_list_df = df[
        (df['Imputed.FDR'] < violin_fdr_threshold) & 
        ((df['log2FC'] < -(violin_log2fc_threshold)) | 
         (df['log2FC'] > violin_log2fc_threshold))
    ]
    pt_level_for_violin = vis.preprocess_for_violin_plot()

    # Customize Plot Section
    with st.expander("Customize Violin Plot"):
        title = st.text_input("Plot Title", st.session_state.violin_config["title"])
        label_size = st.slider(
            "Label Font Size", 
            10, 30, 
            st.session_state.violin_config["label_size"]
        )
        
        # Dynamic color pickers for group-based coloring
        # Use Plotly's qualitative palette for default colors
        group_colors = {}
        unique_groups = vis.annotation_info['Group'].unique()
        default_colors = px.colors.qualitative.Plotly[:len(unique_groups)]

        for i, group in enumerate(unique_groups):
            # Assign a unique default color to each group from Plotly's palette
            default_color = default_colors[i] if i < len(default_colors) else f"#%06x" % (0xFFFFFF & hash(group))
            group_colors[group] = st.color_picker(
                f"Color for {group}",
                st.session_state.violin_config["colors"].get(group, default_color)
            )
        
        # Buttons for Apply and Reset
        apply_col, reset_col = st.columns([1, 1])
        with apply_col:
            if st.button("Apply Changes", key='violin_apply_changes'):
                st.session_state.violin_config.update({
                    "title": title,
                    "label_size": label_size,
                    "colors": group_colors
                })
                st.toast("Plot customization updated!", icon="✅")
        with reset_col:
            if st.button("Reset to Defaults", key='violin_reset_defaults'):
                st.session_state.violin_config = DEFAULT_VIOLIN_CONFIG.copy()
                st.toast("Reset to default settings!", icon="🔄")

    # Protein Selection
    custom_row_select = st.checkbox('Choose custom entries', key='violin_custom_select')
    if custom_row_select:
        subset_df = df[['Protein', 'log2FC', 'Imputed.FDR', 'Gene Name', 'Protein Description']]
        selection = dataframe_with_selections(subset_df, "violin_custom_df_select")
        with st.expander("Your selection"):
            st.write(selection)
        selected_df = pt_level_for_violin[pt_level_for_violin['Protein'].isin(list(selection['Protein']))]
        if len(selected_df) <= 1:
            st.info("Please select at least two inputs to continue.")
        else:
            fig = vis.plot_violin_with_subplots(
                selected_df['Protein'].to_list(), 
                st.session_state.violin_config
            )
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["violin_plot"] = fig
    else:
        dep_list_df['abs_log2FC'] = dep_list_df['log2FC'].abs()
        top10_by_fc = dep_list_df.sort_values(by='abs_log2FC', ascending=False).head(10)
        top_10_df = pt_level_for_violin[pt_level_for_violin['Protein'].isin(list(top10_by_fc['Protein'])[:10])]

        with st.expander("Top 10 DEP"):
            st.dataframe(top_10_df)
        if len(top_10_df) == 0:
            st.info("Your selection has resulted in an empty dataframe. Please adjust thresholds and try again.")
        else:
            fig = vis.plot_violin_with_subplots(
                top_10_df['Protein'].to_list(), 
                st.session_state.violin_config
            )
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["violin_plot"] = fig
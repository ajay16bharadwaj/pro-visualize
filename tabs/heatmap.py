from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_HEATMAP_CONFIG
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

@safe_tab_execution("Heatmap")
@validate_inputs(analysis_status=True, protein_status=True, annotation_status=True)
def render_heatmap(vis, figures_dict, analysis_status, protein_status, annotation_status):
    """Render Heatmap Tab and update figures_dict."""

    # Initialize session state for heatmap_config
    if "heatmap_config" not in st.session_state:
        st.session_state.heatmap_config = DEFAULT_HEATMAP_CONFIG.copy()

    # Preprocess data for heatmap
    heatmap_info = vis.dep_info.copy()  # Differential expression info
    heatmap_condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL)  # Comparisons
    pt_level = vis.preprocess_for_heatmaps()  # Protein-level data for plotting

    # UI for threshold inputs
    col1, col2 = st.columns(2)
    with col1:
        heatmap_log2fc_threshold = st.slider(
            'Log2 Fold Change Threshold',
            0.0, heatmap_info['log2FC'].max(),
            value=0.6, key='heatmap_fc'
        )
    with col2:
        heatmap_fdr_threshold = st.slider(
            'Imputed FDR Threshold',
            0.0, 1.0,
            value=0.05, key='heatmap_fdr'
        )

    # Comparison selection
    heatmap_comparison_input = st.selectbox(
        'Which Comparison',
        heatmap_condition_groups,
        index=0, key="heatmap_comparison_input"
    )
    df = heatmap_condition_groups[heatmap_comparison_input]

    # Filter top DEP proteins based on thresholds
    dep_list_df = df[
        (df['Imputed.FDR'] < heatmap_fdr_threshold) &
        ((df['log2FC'] < -heatmap_log2fc_threshold) | (df['log2FC'] > heatmap_log2fc_threshold))
    ]

    # Display DEP list
    with st.expander("Differentially expressed proteins"):
        st.dataframe(dep_list_df)

    # Custom selection toggle
    custom_row_select = st.checkbox('Choose custom entries', key='heatmap_custom_select')

        # Heatmap customization section
    with st.expander("Customize Heatmap"):
        title = st.text_input("Plot Title", st.session_state.heatmap_config["title"])
        x_label = st.text_input("X-axis Label", st.session_state.heatmap_config["x_label"])
        y_label = st.text_input("Y-axis Label", st.session_state.heatmap_config["y_label"])
        color_gradient = st.selectbox(
            "Color Gradient",
            ["viridis", "plasma", "cividis", "coolwarm", "inferno"],
            index=["viridis", "plasma", "cividis", "coolwarm", "inferno"].index(st.session_state.heatmap_config["color_gradient"]),
            key='heatmap_tab_color_gradient_select'
        )
        height = st.slider("Plot Height", 400, 1000, st.session_state.heatmap_config["height"], key='heatmap_tab_plot_height_slider')
        width = st.slider("Plot Width", 400, 1000, st.session_state.heatmap_config["width"], key='heatmap_tab_plot_width_slider')
        title_font_size = st.slider(
            "Title Font Size",
            10, 40,
            st.session_state.heatmap_config["title_font_size"],
            key='heatmap_tab_title_font_slider'
        )
        x_label_font_size = st.slider(
            "X-axis Font Size",
            10, 30,
            st.session_state.heatmap_config["x_label_font_size"],
            key='heatmap_tab_x_label_font_slider'
        )
        y_label_font_size = st.slider(
            "Y-axis Font Size",
            10, 30,
            st.session_state.heatmap_config["y_label_font_size"],
            key='heatmap_tab_y_label_font_slider'
        )

        # Update heatmap config

    # Buttons for applying changes and resetting to defaults
    apply_col, reset_col = st.columns([1, 1])

    with apply_col:
        if st.button("Apply Changes", key='heatmap_config_apply'):
            st.session_state.heatmap_config.update({
                "title": title,
                "x_label": x_label,
                "y_label": y_label,
                "color_gradient": color_gradient,
                "height": height,
                "width": width,
                "title_font_size": title_font_size,
                "x_label_font_size": x_label_font_size,
                "y_label_font_size": y_label_font_size
            })
            st.toast("Heatmap customization updated!", icon="✅")

    with reset_col:
        if st.button("Reset to Defaults", key='heatmap_config_reset_default'):
            st.session_state.heatmap_config = DEFAULT_HEATMAP_CONFIG.copy()
            st.toast("Reset to default settings!", icon="🔄")

    if custom_row_select:
        # Custom protein selection
        subset_df = df[['Protein', 'log2FC', 'Imputed.FDR', 'Gene Name', 'Protein Description']]
        selection = dataframe_with_selections(subset_df, "heatmap_custom_df_select")
        with st.expander("Your selection"):
            st.write(selection)

        selected_df = pt_level[pt_level['Protein'].isin(list(selection['Protein']))]
        if len(selected_df) <= 1:
            st.info("Please select at least two inputs to continue.")
        else:
            fig = vis.plot_heatmap(
                selected_df,
                st.session_state.heatmap_config
            )
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["heatmap"] = fig
    else:
        # Top 10 DEP proteins based on absolute log2FC
        dep_list_df['abs_log2FC'] = dep_list_df['log2FC'].abs()
        top10_by_fc = dep_list_df.sort_values(by='abs_log2FC', ascending=False).head(10)
        top_10_df = pt_level[pt_level['Protein'].isin(list(top10_by_fc['Protein'])[:10])]

        with st.expander("Top 10 DEP"):
            st.dataframe(top_10_df)

        if len(top_10_df) == 0:
            st.info("Your selection has resulted in an empty dataframe. Please select different thresholds and try again.")
        else:
            fig = vis.plot_heatmap(
                top_10_df,
                st.session_state.heatmap_config
            )
            st.plotly_chart(fig, theme="streamlit", use_container_width=True)
            figures_dict["heatmap"] = fig


from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_PCA_BY_ANNOTATION_CONFIG, DEFAULT_UMAP_CONFIG, DEFAULT_TSNE_CONFIG  
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

#clustering plots - Tab Code. Need documentation for what the plot is? 
#@safe_tab_execution("Clustering")
@validate_inputs(protein_status=True, annotation_status=True)
def render_clustering(vis, figures_dict, protein_status, annotation_status, global_config, **kwargs):
    """Render clustering plots Tab and update figures_dict."""
    clust_tab1, clust_tab2, clust_tab3 = st.tabs(["PCA", "UMAP", "T-SNE"])
    with clust_tab1: 
        st.write(' This tab will be used for Clustering')
        # Placeholder for four example plots

        #checks

        #PCA computation. 
        vis.preprocess_for_pca()
        #st.write(vis.protein_data_for_pca)
        # Initialize session state for PCA config
        if "pca_config" not in st.session_state:
            st.session_state.pca_config = DEFAULT_PCA_BY_ANNOTATION_CONFIG.copy()

        # Create a temporary config for user inputs
        user_inputs = st.session_state.pca_config.copy()

        # User inputs for customizing the PCA plot
        with st.expander("Customize PCA Plot"):
            user_inputs["title"] = st.text_input(
                "Plot Title", 
                user_inputs["title"], 
                key="pca_plot_title"
            )
            user_inputs["x_label"] = st.text_input(
                "X-axis Label", 
                user_inputs["x_label"], 
                key="pca_x_label"
            )
            user_inputs["y_label"] = st.text_input(
                "Y-axis Label", 
                user_inputs["y_label"], 
                key="pca_y_label"
            )
            user_inputs["marker_size"] = st.slider(
                "Marker Size", 
                5, 20, 
                user_inputs["marker_size"], 
                key="pca_marker_size"
            )
            user_inputs["marker_symbol"] = st.selectbox(
                "Marker Symbol", 
                ["circle", "star", "cross", "diamond", "square", "plus"], 
                index=["circle", "star", "cross", "diamond", "square", "plus"].index(user_inputs.get("marker_symbol", "circle")
                ),                                                                                     
                key="pca_marker_symbol"
            )
            user_inputs["jitter"] = st.slider(
                "Jitter (Overlap Reduction)", 
                0.0, 1.0, 
                user_inputs["jitter"], 
                key="pca_jitter"
            )
            user_inputs["width"] = st.slider(
                "Plot Width", 
                400, 1200, 
                user_inputs["width"], 
                key="pca_width"
            )
            user_inputs["height"] = st.slider(
                "Plot Height", 
                400, 1200, 
                user_inputs["height"], 
                key="pca_height"
            )

        # Apply and Reset buttons
        apply_col, reset_col = st.columns([1, 1])
        with apply_col:
            if st.button("Apply Changes", key="pca_apply_changes"):
                st.session_state.pca_config = user_inputs.copy()
                st.toast("PCA plot configuration updated!", icon="✅")
        with reset_col:
            if st.button("Reset to Default", key="pca_reset_default"):
                st.session_state.pca_config = DEFAULT_PCA_BY_ANNOTATION_CONFIG.copy()
                st.toast("PCA plot configuration reset to defaults!", icon="🔄")

        st.write(config=st.session_state.pca_config)
        
        # Generate the PCA plot using updated configuration and global group colors
        pca_plot_by_annotation = vis.plot_pca_by_annotation(
            config=st.session_state.pca_config,
            group_column="Group",
            color_discrete_map=global_config["group_colors"]
        )
        figures_dict["pca_plot_by_annotation"] = pca_plot_by_annotation

        #pca_plot_by_annotation = vis.plot_pca_by_annotation()
        pca_plot_by_clusters = vis.plot_pca_with_clusters_plotly()
        hierarchial_clustering_dendogram = vis.plot_vertical_dendrogram()
        cluster_assignment_table = vis.create_cluster_assignment_table()


    #arranging the page display
        col1, col2 = st.columns(2) 
                
        with col1:
            st.plotly_chart(pca_plot_by_annotation, use_container_width=True)  # Plot 1 in the first column
            st.write("sample")
        with col2:
            st.plotly_chart(pca_plot_by_clusters, use_container_width=True)  # Plot 2 in the second column

        # # Row 2: Another two plots side by side
        col3, col4 = st.columns(2)
        
        with col3:
            st.pyplot(hierarchial_clustering_dendogram, use_container_width=True)  # Plot 3 in the first column
        with col4:
            st.write("Cluster Assignment Table")
            st.dataframe(cluster_assignment_table, use_container_width=True)

    with clust_tab2:
        #st.write("will integrate umap here")
        st.write("UMAP Plot")

        # Initialize session state for UMAP config
        if "umap_config" not in st.session_state:
            st.session_state.umap_config = DEFAULT_UMAP_CONFIG.copy()
            st.session_state.umap_config["colors"] = global_config["group_colors"].copy()

        # Synchronize session state colors with global config
        for group, color in global_config["group_colors"].items():
            if group not in st.session_state.umap_config["colors"]:
                st.session_state.umap_config["colors"][group] = color
            elif st.session_state.umap_config["colors"][group] != color:
                st.session_state.umap_config["colors"][group] = color

        # Temporary UMAP config for user input
        temp_config = st.session_state.umap_config.copy()

        # Customization block
        with st.expander("Customize UMAP Plot"):
            temp_config["title"] = st.text_input(
                "Plot Title", temp_config["title"], key="umap_title"
            )
            temp_config["marker_size"] = st.slider(
                "Marker Size", 5, 20, temp_config["marker_size"], key="umap_marker_size"
            )
            temp_config["width"] = st.slider(
                "Plot Width", 400, 1200, temp_config["width"], key="umap_width"
            )
            temp_config["height"] = st.slider(
                "Plot Height", 400, 1200, temp_config["height"], key="umap_height"
            )

            # Apply and Reset buttons
            apply_col, reset_col = st.columns([1, 1])
            with apply_col:
                if st.button("Apply Changes", key="umap_apply_changes"):
                    st.session_state.umap_config = temp_config.copy()
                    st.toast("UMAP plot configuration updated!", icon="✅")
            with reset_col:
                if st.button("Reset to Defaults", key="umap_reset_default"):
                    st.session_state.umap_config = DEFAULT_UMAP_CONFIG.copy()
                    st.session_state.umap_config["colors"] = global_config["group_colors"].copy()
                    st.toast("UMAP plot configuration reset to defaults!", icon="🔄")

        # Generate UMAP plot
        umap_plot = vis.plot_umap(config=st.session_state.umap_config)
        st.plotly_chart(umap_plot, use_container_width=True)

    with clust_tab3:
        st.write("t-SNE Plot")

        # Initialize session state for t-SNE config
        if "tsne_config" not in st.session_state:
            st.session_state.tsne_config = DEFAULT_TSNE_CONFIG.copy()
            st.session_state.tsne_config["colors"] = global_config["group_colors"].copy()

        # Synchronize session state colors with global config
        for group, color in global_config["group_colors"].items():
            if group not in st.session_state.tsne_config["colors"]:
                st.session_state.tsne_config["colors"][group] = color
            elif st.session_state.tsne_config["colors"][group] != color:
                st.session_state.tsne_config["colors"][group] = color

        # Temporary t-SNE config for user input
        temp_config = st.session_state.tsne_config.copy()

        # Customization block
        with st.expander("Customize t-SNE Plot"):
            temp_config["title"] = st.text_input(
                "Plot Title", temp_config["title"], key="tsne_title"
            )
            temp_config["marker_size"] = st.slider(
                "Marker Size", 5, 20, temp_config["marker_size"], key="tsne_marker_size"
            )
            temp_config["width"] = st.slider(
                "Plot Width", 400, 1200, temp_config["width"], key="tsne_width"
            )
            temp_config["height"] = st.slider(
                "Plot Height", 400, 1200, temp_config["height"], key="tsne_height"
            )
            temp_config["show_labels"] = st.checkbox(
                "Show Sample Labels", temp_config["show_labels"], key="tsne_show_labels"
            )

            # Apply and Reset buttons
            apply_col, reset_col = st.columns([1, 1])
            with apply_col:
                if st.button("Apply Changes", key="tsne_apply_changes"):
                    st.session_state.tsne_config = temp_config.copy()
                    st.toast("t-SNE plot configuration updated!", icon="✅")
            with reset_col:
                if st.button("Reset to Defaults", key="tsne_reset_default"):
                    st.session_state.tsne_config = DEFAULT_TSNE_CONFIG.copy()
                    st.session_state.tsne_config["colors"] = global_config["group_colors"].copy()
                    st.toast("t-SNE plot configuration reset to defaults!", icon="🔄")

        # Generate t-SNE plot
        tsne_plot = vis.plot_tsne(config=st.session_state.tsne_config)
        st.plotly_chart(tsne_plot, use_container_width=True)
    
    
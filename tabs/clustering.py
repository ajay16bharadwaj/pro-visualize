from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd
from utils.helpers import dataframe_with_selections

#clustering plots - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("Clustering")
@validate_inputs(protein_status=True, annotation_status=True)
def render_clustering(vis, figures_dict, protein_status, annotation_status, **kwargs):
    """Render clustering plots Tab and update figures_dict."""
    clust_tab1, clust_tab2, clust_tab3 = st.tabs(["PCA", "UMAP", "T-SNE"])
    with clust_tab1: 
        st.write(' This tab will be used for Clustering')
        # Placeholder for four example plots

        #checks

        #PCA computation. 
        vis.preprocess_for_pca()
        #st.write(vis.protein_data_for_pca)
        pca_plot_by_annotation = vis.plot_pca_by_annotation()
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
        umap_plot = vis.plot_umap()
        st.plotly_chart(umap_plot, use_container_width=True)

    with clust_tab3:
        tsne_plot = vis.plot_tsne()
        st.plotly_chart(tsne_plot, use_container_width=True)
    
    
from matplotlib import pyplot as plt
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from compute_functions import *
import numpy as np
import dash_bio
from sklearn.preprocessing import MinMaxScaler
from visualization import ProteinVisualization

##################
#CONSTANTS
#################
COL_DEPLABEL = 'Label2'
COL_DEPSIGNIF = 'Imputed.FDR'
color_map = {
        'Down-regulated': '#0072B2',
        'Significant, no change': '#F0E442',
        'Up-regulated': '#D55E00',
        'Non-significant': '#999999'
    }

######################################
# functions
######################################
@st.cache_data
def load_data(path):
    df = pd.read_csv(path,sep = '\t')
    return df



## should I move this to compute or helper functions? 
def dataframe_with_selections(df, custom_key_name):
    df_with_selections = df.copy()
    df_with_selections.insert(0, "Select", False)

    # Get dataframe row-selections from user with st.data_editor
    edited_df = st.data_editor(
        df_with_selections,
        hide_index=True,
        column_config={"Select": st.column_config.CheckboxColumn(required=True)},
        disabled=df.columns,
        key=custom_key_name
    )

    # Filter the dataframe using the temporary column, then drop the column
    selected_rows = edited_df[edited_df.Select]
    return selected_rows.drop('Select', axis=1)

#######################################
# PAGE SETUP
#######################################

st.set_page_config(page_title="ProEpic Report", page_icon=":bar_chart:", layout="wide")

st.title("Pro-Visualize")
st.markdown("_Prototype v0.4.1_")

analysis_status = False
protein_level_status = False
annotation_status = False

with st.sidebar:

    st.header("Inputs")
    analysis_upload = st.file_uploader("Analysis Input")

    if analysis_upload is None:
        st.info(" Upload Differentially expressed file through config", icon="ℹ️")

    protein_level_upload = st.file_uploader("Protein Level Input")

    if protein_level_upload is None:
        st.info( " Upload the Protein Level input through config", icon="ℹ️")

    annotation_file_upload = st.file_uploader("Upload your metadata/Annotation file")

    if annotation_file_upload is None:
        st.info( " Upload the metadata/Annotation file through config", icon="ℹ️")
        


#######################################
# DATA LOADING
#######################################

# Initialize dictionary to store the figures
figures_dict = {}
#intialize class 
vis = ProteinVisualization() 

with st.container():

    col1, col2, col3 = st.columns(3)

    with col1:
        with st.expander("Analysis Input Preview"):
            if analysis_upload is not None:
                #analysis_df = load_data(analysis_upload) # type: ignore
                vis.load_dep_data(analysis_upload) # type: ignore
                #vis.dep_info = analysis_df.copy()
                #vis.analysis_with_annotation = analysis_df.copy()
                st.dataframe(vis.dep_info)
                analysis_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")
    with col2:
        with st.expander("Annotation Preview"):
            if annotation_file_upload is not None:
                #annotation_file_upload = load_data(annotation_file_upload)
                vis.load_annotation(annotation_file_upload, "Level3", "attribute_ExperimentalGroup" ) # type: ignore #parameter values later on. 
                st.dataframe(vis.annotation_info)
                annotation_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")

    with col3:
        with st.expander("Protein Input Preview"):
            if protein_level_upload is not None:
                #protein_level_upload = load_data(protein_level_upload)
                vis.load_protein_data(protein_level_upload) # type: ignore
                #upload pt level with uniport annotation - to save time. 
                vis.protein_data_annotated = vis.protein_data.copy()
                st.dataframe(vis.protein_data_annotated)
                protein_level_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")

    #  if analysis_status and protein_level_status and annotation_status:
    #     with st.spinner("Fetching UniProt information..."):
    #             vis.get_uniprot_info() 
    #             vis.merge_uniprot2proteome()  
    #             st.success("UniProt information fetched successfully!")
 

#initializing tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["Volcano Plot", "Heat Maps", "Violin Plot", "Quantification", "Clustering"])

#volcano plot - to set interactive FDR and FC thresholds
with tab1:
    if analysis_status and annotation_file_upload:
        #st.write('this tab will be used for volcano plot')
        volcano_info = vis.dep_info.copy() # type: ignore
        condition_groups = vis.volcano_preprocess(COL_DEPLABEL) # type: ignore
        volcano_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, volcano_info['log2FC'].max(), value=0.6, key='volcano_fc') # type: ignore
        volcano_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='volcano_fdr')
        volcano_comparison_input = st.selectbox('Which comparison', condition_groups, index=1, key="volcano_comparison_input" )   
        #st.dataframe(condition_groups[volcano_comparison_input])     
        fig = vis.plot_volcano(condition_groups[volcano_comparison_input], volcano_log2fc_threshold, volcano_fdr_threshold)
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)
        figures_dict["Volcano Plot"] = fig
        

#heatmaps
#default to top 10 DEP, and options FDR and FC values. Give the option to select custom proteins as well. 
with tab2:
    if protein_level_status and analysis_status:
        st.write('this tab will be used for the heatmaps')
        if protein_level_status and analysis_status and annotation_status:

            pt_level = vis.preprocess_for_heatmaps()
            heatmap_info = vis.dep_info.copy() #type:ignore
            heatmap_condition_groups = vis.volcano_preprocess(COL_DEPLABEL) # type: ignore
            heatmap_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, heatmap_info['log2FC'].max(), value=0.6, key='heatmap_fc') # type: ignore
            heatmap_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='heatmap_fdr')
            heatmap_comparison_input = st.selectbox('Which comparison', heatmap_condition_groups, index=1, key="heatmap_comparison_input" )

            #getting the top 10 differentially expressed proteins and plotting that 
            df = heatmap_condition_groups[heatmap_comparison_input]
            dep_list_df = df[(df['Imputed.FDR'] < heatmap_fdr_threshold) & ((df['log2FC'] < heatmap_log2fc_threshold) | (df['log2FC'] > heatmap_log2fc_threshold)) ]

            #protein level preprocessing for violin plot
            pt_level_for_heatmap = vis.preprocess_for_heatmaps()
            
            with st.expander("Differentially expressed proteins"):
                st.dataframe(dep_list_df)

            custom_row_select = st.checkbox(' Choose custom entries ', key='heatmap_custom_select')

            if custom_row_select:
                selection = dataframe_with_selections(df, "heatmap_custom_df_select")
                with st.expander("Your selection"):
                    st.write(selection)

                
                selected_df = pt_level_for_heatmap[pt_level_for_heatmap['Protein'].isin(list(selection['Protein']))]
                if (len(selected_df) <= 1):
                    st.info("Please select at-least two inputs to continue")
                else:
                    fig = vis.plot_heatmap(selected_df)
                    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
                    figures_dict["Heatmap"] = fig

            else:
                dep_list_df['abs_log2FC'] = dep_list_df['log2FC'].abs()  # Create a new column with the absolute values of log2FC
                top10_by_fc = dep_list_df.sort_values(by='abs_log2FC', ascending=False).head(10)
                top_10_df = pt_level_for_heatmap[pt_level_for_heatmap['Protein'].isin(list(top10_by_fc['Protein'])[:10])]

                with st.expander("Top 10 DEP"):
                    st.dataframe(top_10_df)

                if (len(top_10_df) == 0):
                    st.info("Your selection has resulted in empty dataframe, please select different thresholds and try again")
                else:
                    fig = vis.plot_heatmap(top_10_df)
                    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
                    figures_dict["Heatmap"] = fig

#this tab is for the violin plots. 
with tab3:
    st.write('this tab will be used for violin plot')
    if protein_level_status and analysis_status and annotation_status:
        violin_info = vis.dep_info.copy() #type:ignore
        violin_condition_groups = vis.volcano_preprocess(COL_DEPLABEL) # type: ignore
        violin_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, violin_info['log2FC'].max(), value=0.6, key='violin_fc') # type: ignore
        violin_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='violin_fdr')
        violin_comparison_input = st.selectbox('Which comparison', violin_condition_groups, index=1, key="violin_comparison_input" )

        #getting the top 10 differentially expressed proteins and plotting that 
        df = violin_condition_groups[violin_comparison_input]
        dep_list_df = df[(df['Imputed.FDR'] < violin_fdr_threshold) & ((df['log2FC'] < violin_log2fc_threshold) | (df['log2FC'] > violin_log2fc_threshold)) ]

        #protein level preprocessing for violin plot
        pt_level_for_violin = vis.preprocess_for_violin_plot()
        
        with st.expander("Differentially expressed proteins"):
            st.dataframe(dep_list_df)

        custom_row_select = st.checkbox(' Choose custom entries ', key='violin_custom_select')

        if custom_row_select:
            selection = dataframe_with_selections(df, "violin_custom_df_select")
            with st.expander("Your selection"):
                st.write(selection)

            
            selected_df = pt_level_for_violin[pt_level_for_violin['Protein'].isin(list(selection['Protein']))]
            if (len(selected_df) <= 1):
                st.info("Please select at-least two inputs to continue")
            else:
                fig = vis.plot_violin_with_subplots(selected_df['Protein'].to_list())
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)

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


#this tab is used for quantification. 
with tab4:
    st.write('this tab will be used for quantification plots')
    if protein_level_status and annotation_status:
        quant_tab1, quant_tab2, quant_tab3, quant_tab4, quant_tab5 = st.tabs(['Protein Per Sample', 'Protein Overlap', 'Protein Intensity Density', 'Correlation Matrix', 'Protein Rank Order'])
        
        with quant_tab1: 
            plot_proteins_per_sample =  vis.plot_proteins_per_sample()
            st.plotly_chart(plot_proteins_per_sample)

        with quant_tab2: 
            plot_protein_overlap = vis.plot_protein_overlap()
            st.plotly_chart(plot_protein_overlap)

        with quant_tab3:
            plot_intensity_density = vis.plot_intensity_density()
            st.plotly_chart(plot_intensity_density)

        with quant_tab4:
            plot_correlation_matrix = vis.plot_correlation_matrix()
            st.plotly_chart(plot_correlation_matrix)

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

            else:
                protein_rank_order_plot = vis.plot_protein_rank_order()
                st.plotly_chart(protein_rank_order_plot, use_container_width=True)


        


#this tab will be used for PCA 
with tab5:

    if protein_level_status and annotation_status:
        clust_tab1, clust_tab2, clust_tab3 = st.tabs(["PCA", "UMAP", "T-SNE"])
        with clust_tab1: 
            st.write(' This tab will be used for Clustering')
            # Placeholder for four example plots

            #checks
            if protein_level_status and annotation_status:

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

with st.sidebar:
    create_download_button(figures_dict)

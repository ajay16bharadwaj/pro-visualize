import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from compute_functions import *
import numpy as np
import dash_bio
from sklearn.preprocessing import MinMaxScaler

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
def load_data(path: str):
    df = pd.read_csv(path,sep = '\t')
    return df


#function to plot the volcano plot? 
def plot_volcano(df, FC, p_val):
    # Creating the plot
    df['sig'] = -np.log10(df[COL_DEPSIGNIF])

    df['label'] = 'Non-significant'
    df.loc[df[COL_DEPSIGNIF] <= p_val, 'label'] = 'Significant, no change'
    df.loc[(df[COL_DEPSIGNIF] <= p_val) & (df['log2FC'] > FC), 'label'] = 'Up-regulated'
    df.loc[(df[COL_DEPSIGNIF] <= p_val) & (df['log2FC'] < -FC), 'label'] = 'Down-regulated'

    fig = px.scatter(df, x='log2FC', y='sig', color='label', color_discrete_map=color_map, hover_data=['Gene Name', 'Protein Description'])
        
        # Customizing the layout
    fig.update_layout(title=f'Volcano Plot', xaxis_title="log2 fold change", yaxis_title="-log10(FDR)", legend_title_text='Category')
        
        # Adding horizontal and vertical lines for cutoffs
    fig.add_shape(type="line", x0=-FC, y0=0, x1=-FC, y1=df['sig'].max(), line=dict(color="RoyalBlue", width=2, dash="dot"))
    fig.add_shape(type="line", x0=FC, y0=0, x1=FC, y1=df['sig'].max(), line=dict(color="RoyalBlue", width=2, dash="dot"))
    fig.add_shape(type="line", x0=df['log2FC'].min(), y0=-np.log10(p_val), x1=df['log2FC'].max(), y1=-np.log10(p_val),line=dict(color="Red", width=2, dash="dot"))

    return fig    


def plot_heatmap(df):

    #need to set the index and scale it before plotting. 
    df = df.set_index('Protein')
    #df = df.drop('Protein', axis=1)

    scaler = MinMaxScaler()
    df1_norm = pd.DataFrame(scaler.fit_transform(df), columns=df.columns, index=df.index)


    fig = dash_bio.Clustergram(
    data=df1_norm,
    column_labels=list(df1_norm.columns.values),
    row_labels=list(df1_norm.index),
    height=600,
    width=600
    )

    return fig

def dataframe_with_selections(df):
    df_with_selections = df.copy()
    df_with_selections.insert(0, "Select", False)

    # Get dataframe row-selections from user with st.data_editor
    edited_df = st.data_editor(
        df_with_selections,
        hide_index=True,
        column_config={"Select": st.column_config.CheckboxColumn(required=True)},
        disabled=df.columns,
    )

    # Filter the dataframe using the temporary column, then drop the column
    selected_rows = edited_df[edited_df.Select]
    return selected_rows.drop('Select', axis=1)

#######################################
# PAGE SETUP
#######################################

st.set_page_config(page_title="ProEpic Report", page_icon=":bar_chart:", layout="wide")

st.title("Auto Report Generator")
st.markdown("_Prototype v0.4.1_")

analysis_status = False
protein_level_status = False

with st.sidebar:

    st.header("Inputs")
    analysis_upload = st.file_uploader("Analysis Input")

    if analysis_upload is None:
        st.info(" Upload Differentially expressed file through config", icon="ℹ️")

    protein_level_upload = st.file_uploader("Protein Level Input")

    if protein_level_upload is None:
        st.info( " Upload the Protein Level input through config", icon="ℹ️")
        


#######################################
# DATA LOADING
#######################################

with st.container():

    col1, col2 = st.columns(2)

    with col1:
        with st.expander("Analysis Input Preview"):
            if analysis_upload is not None:
                analysis_df = load_data(analysis_upload) # type: ignore
                st.dataframe(analysis_df)
                analysis_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")

    with col2:
        with st.expander("Protein Input Preview"):
            if protein_level_upload is not None:
                protein_level_upload = load_data(protein_level_upload) # type: ignore
                st.dataframe(protein_level_upload)
                protein_level_status = True
            else:
                st.info(" Preview available after file upload", icon="ℹ️")

#initializing tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["Volcano Plot", "Heat Maps", "Violin Plot", "Quantification", "LLM"])

#volcano plot - to set interactive FDR and FC thresholds
with tab1:
    if analysis_status:
        volcano_info = analysis_df.copy()        
        #don't need this anymore because I'm uploading file with uniprot information
        #uniprot_annotation = get_uniprot_info(volcano_info)
        #volcano_info = merge_uniprot2proteome(volcano_info, uniprot_annotation)
        #print(volcano_info.head())
        condition_level = sorted(volcano_info[COL_DEPLABEL].unique())
        volcano_info[COL_DEPLABEL] = pd.Categorical(volcano_info[COL_DEPLABEL], categories=condition_level, ordered=True)
        condition_groups = {condition: group for condition, group in volcano_info.groupby(COL_DEPLABEL)}
        log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, volcano_info['log2FC'].max(), value=0.6, key='volcano_fc')
        fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='volcano_fdr')
        comparison_input = st.selectbox('Which comparison', condition_groups, index=1 )
        fig = plot_volcano(condition_groups[comparison_input], log2fc_threshold, fdr_threshold) # type: ignore
        st.plotly_chart(fig, theme="streamlit", use_container_width=True)

#heatmaps
#default to top 10 DEP, and options FDR and FC values. Give the option to select custom proteins as well. 
with tab2:
    if protein_level_status and analysis_status:

        pt_level = protein_level_upload.copy() # type: ignore

        #cleaning pt_level
        # Replace Inf, -Inf with NaN first (if you specifically only want to replace Inf values)
        pt_level.replace([np.inf, -np.inf], np.nan, inplace=True)

        # Then replace all NaN values with 0
        pt_level.fillna(0, inplace=True)

        volcano_info = analysis_df.copy()
        condition_level = sorted(volcano_info[COL_DEPLABEL].unique())
        volcano_info[COL_DEPLABEL] = pd.Categorical(volcano_info[COL_DEPLABEL], categories=condition_level, ordered=True)
        condition_groups = {condition: group for condition, group in volcano_info.groupby(COL_DEPLABEL)}

        heatmap_comparison_input = st.selectbox(' Which comparison', condition_groups, index=1 )

        heatmap_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, volcano_info['log2FC'].max(), value=0.6, key='heatmap_fc')
        heatmap_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='heatmap_fdr')

        #Need to get top ten DEP's and plot that 
        df = condition_groups[heatmap_comparison_input] # type: ignore
        dep_list_df = df[(df['Imputed.FDR'] < heatmap_fdr_threshold) & ((df['log2FC'] < heatmap_log2fc_threshold) | (df['log2FC'] > heatmap_log2fc_threshold)) ]
        
        with st.expander("Differentially expressed proteins"):
            st.dataframe(dep_list_df)

        custom_row_select = st.checkbox(' Choose custom entries ')

        if custom_row_select:
            selection = dataframe_with_selections(df)
            with st.expander("Your selection"):
                st.write(selection)

            selected_df = pt_level[pt_level['Protein'].isin(list(selection['Protein']))]
            if (len(selected_df) <= 1):
                st.info("Please select at-least two inputs to continue")
            else:
                fig = plot_heatmap(selected_df)
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)

        else:
            dep_list_df['abs_log2FC'] = dep_list_df['log2FC'].abs()  # Create a new column with the absolute values of log2FC
            top10_by_fc = dep_list_df.sort_values(by='abs_log2FC', ascending=False).head(10)
            top_10_df = pt_level[pt_level['Protein'].isin(list(top10_by_fc['Protein'])[:10])]

            with st.expander("Top 10 DEP"):
                st.dataframe(top_10_df)

            if (len(top_10_df) == 0):
                st.info("Your selection has resulted in empty dataframe, please select different thresholds and try again")
            else:
                fig = plot_heatmap(top_10_df)
                st.plotly_chart(fig, theme="streamlit", use_container_width=True)

#this tab is for the violin plots. 
with tab3:
    st.write('this tab will be used for violin plot')


with tab5:
    st.write(' This tab will be used for an LLM')

    
from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD 
import streamlit as st
import pandas as pd
from io import BytesIO
from utils.helpers import dataframe_with_selections

#uniprot_annotation  - get uniprot annotation with protein level file upload
@safe_tab_execution("uniprot_annotation")
def render_get_uniprot_annotation(vis, **kwargs):
    """Render uniprot annotation """
    st.write(" If your protein data is not Uniprot annotated, please run this first")
    unannoteted_protein_level_upload = st.file_uploader("Unannotated Protein Level Input")

    if unannoteted_protein_level_upload is not None:
        #protein_level_upload = load_data(protein_level_upload)
        vis.load_protein_data(unannoteted_protein_level_upload) # type: ignore

        # Display a spinner with a message while get_uniprot_info is running
        with st.spinner("Annotating proteins... Please wait."):
            vis.get_uniprot_info()
            annotated_protein_data = vis.merge_uniprot2proteome(vis.protein_data)

        st.dataframe(annotated_protein_data)

        protein_annotated_data_download = annotated_protein_data.to_csv(sep='\t', index=False)

        st.download_button(
            label="Download data as TSV",
            data=protein_annotated_data_download,
            file_name='protein_level_annotated.txt',
            mime='text/csv',
            help="Click to download the current data as a TSV/TXT file"
        )
    
    
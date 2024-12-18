import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import time
from io import StringIO
from requests.adapters import HTTPAdapter, Retry
import re
import streamlit as st
import plotly.io as pio
from io import BytesIO
import base64



def get_uniprot_entries(protein_names, split_list, cols):

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    entries = []

    for protein_list in split_list:
        query_proteins = "+OR+".join(protein_list)
        url = f'https://rest.uniprot.org/uniprotkb/search?fields={cols}&size=300&format=tsv&query={query_proteins}'

        for batch,total in get_batch(url):
            for line in batch.text.splitlines()[1:]:
                line = line.split("\t")
                if line[0] not in protein_names:
                    continue
                if line not in entries:
                    entries.append(line)

    return entries

def get_uniprot_info(pt_level):

    NoSplitRequired = None
    proteins = list(set(pt_level['Protein'].tolist()))
    protein_names = []

    #removes fasta descriptions if present
    for protein_name in proteins:
        pt = _get_uniprot_ids(protein_name)
        protein_names.append(pt)

    #splitting protein list into separate lists since the Uniprot API only allows queries under a certain size
    chunk_size = 130
    if len(protein_names) > chunk_size:
        split_list = []
        for i in range(0, len(protein_names), chunk_size):
            split_list.append(protein_names[i:i+chunk_size])
    else:
        NoSplitRequired = True
        split_list = [protein_names]

    #defining columns to retrieve from uniprot: https://www.uniprot.org/help/return_fields
    cols = ["accession", "gene_primary","protein_name","cc_tissue_specificity","go_p","go_c","go_f","cc_subcellular_location"]
    cols = ",".join(cols)

    entries_list = get_uniprot_entries(protein_names, split_list, cols)

    if len(entries_list) != 0:
        #change the column headers
        df_cols = ["ProteinIds","Gene Name", "Protein Description", "Tissue Specificity","Gene ontology (biological process)","Gene ontology (cellular component)","Gene ontology (molecular function)","Subcellular Location[CC]"]
        uniprot_df = pd.DataFrame(entries_list, columns=df_cols)

    else:
        print("Couldn't fetch UniProt Annotation at the moment.")

    return uniprot_df

def merge_uniprot2proteome(pt_level, uniprot_df):

    pt_level['ProteinIds'] = pt_level['Protein'].apply(_get_uniprot_ids)
    pt_level = pd.merge(uniprot_df, pt_level, on="ProteinIds", how="right")
    
    return pt_level

## helper functions
def _get_uniprot_ids(row):
        REGEX_UNIPROTID = re.compile("iRT|([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?")
        result = re.search(REGEX_UNIPROTID, row)
        if result is None:
                pt = row
        else:
                pt = result.group(0)
        return pt


# Function to export figures as HTML
def export_figures_to_html(figures_dict):
    """
    Export the figures to an HTML file.
    Args:
        figures_dict (dict): A dictionary containing the figure names and their Plotly objects.
    Returns:
        str: HTML content of the figures.
    """
    html_content = "<html><head><title>Exported Figures</title></head><body>"
    html_content += "<h1>Exported Figures from the Streamlit App</h1>"

    for fig_name, fig in figures_dict.items():
        html_content += f"<h2>{fig_name}</h2>"
        html_content += pio.to_html(fig, include_plotlyjs='cdn', full_html=False) # type: ignore

    html_content += "</body></html>"
    
    return html_content

# Function to generate a downloadable HTML file
def create_download_button(figures_dict, filename="exported_figures.html"):
    """
    Creates a download button for the exported figures.
    Args:
        figures_dict (dict): A dictionary containing the figure names and their Plotly objects.
        filename (str): The filename for the exported HTML file.
    """
    # Generate the HTML content
    html_content = export_figures_to_html(figures_dict)

    # Convert the HTML content to bytes for download
    b64 = base64.b64encode(html_content.encode()).decode()

    # Create a download button
    st.download_button(
        label="Download all figures as HTML",
        data=BytesIO(html_content.encode()),
        file_name=filename,
        mime="text/html"
    )

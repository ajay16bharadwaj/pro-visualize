from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_GENE_NAME_COLUMN 
import streamlit as st
import pandas as pd
from io import BytesIO
from utils.helpers import dataframe_with_selections
from utils.streamlit_caching import cached_get_all_enrichment
import networkx as nx
from pyvis.network import Network

#functional_annotation plots - Tab Code. Need documentation for what the plot is? 
@safe_tab_execution("functional_annotation")
@validate_inputs(analysis_status=True,protein_status=True, annotation_status=True)
def render_functional_annotation(vis, figures_dict, analysis_status, protein_status, annotation_status, **kwargs):
    """Render functional annotation and update figures_dict."""
    df = vis.dep_info.copy() # type: ignore
    anaysis_condition_groups = vis.volcano_preprocess(DEFAULT_COL_DEPLABEL) # type: ignore
    analysis_log2fc_threshold = st.slider('Log2 Fold Change Threshold', 0.0, df['log2FC'].max(), value=0.6, key='analysis_fc') # type: ignore
    analysis_fdr_threshold = st.slider('Imputed FDR Threshold', 0.0, 1.0, 0.05, key='analysis_fdr')
    analysis_comparison_input = st.selectbox('Which comparison', anaysis_condition_groups, index=0, key="analysis_comparison_input" )
    organism_input = st.selectbox(' Choose your organism', list(vis.organism_dict.keys()), index=0, key="organism_input")

    #choosing comparison based filtered df
    filtered_df = anaysis_condition_groups[analysis_comparison_input]
    #getting the list of differentially expressed proteins for that comparison
    dep_list_df = filtered_df[(filtered_df['Imputed.FDR'] < analysis_fdr_threshold) & ((filtered_df['log2FC'] < -(analysis_log2fc_threshold)) | (filtered_df['log2FC'] > analysis_log2fc_threshold)) ]
    if dep_list_df.empty:
        st.warning("No differentially expressed proteins found for the selected comparison. Please select a different comparison.")
        st.stop()  # Stops the execution here, so the user has to select a valid comparison.
    enrichment_df = None
    source_dict = None

    #getting the terms with this set of significant proteins. 
    ea, go_cc, go_bp, go_mf, kegg, Network_Visualization = st.tabs(['Comprehensive Enrichment Analysis','GO Cellular Component Encrichment', 'GO Biological Process', 'GO Molecular Function', 'KEGG Biological Pathways', 'Network Visualization'])

    with ea:
        enrichment_df, source_dict =  cached_get_all_enrichment(vis, list(dep_list_df[DEFAULT_GENE_NAME_COLUMN]), organism_input)
        ea_manhattan_plot = vis.plot_manhattan(enrichment_df, category_name="Comprehensive Enrichment Analysis")
        st.plotly_chart(ea_manhattan_plot)

    with go_cc:
        if not enrichment_df.empty and source_dict is not None:
            if "GO:CC" in list(source_dict.keys()):
                cc_df = source_dict["GO:CC"]

                if not cc_df.empty:
                    # Display the GO:CC dot plot
                    #st.subheader("GO:CC Dot Plot")

                    # Initialize session state for selected terms if not already set
                    if "cc_selected_terms" not in st.session_state:
                        st.session_state.cc_selected_terms = pd.DataFrame()
                    if "cc_selection_state" not in st.session_state:
                        st.session_state.cc_selection_state = []  # Track selected row indices

                    # Render the plot
                    if st.session_state.cc_selected_terms.empty:
                        # Default plot with top 10 GO terms
                        fig_cc, ax_cc = vis.plot_go_dotplot(
                            cc_df.sort_values(by="q_value").head(10),
                            category_name="Cellular Component"
                        )
                        fig_cc.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                        st.pyplot(fig_cc)
                    else:
                        # Filter `cc_df` to include only the selected terms
                        selected_native_ids = st.session_state.cc_selected_terms["native"].tolist()
                        filtered_cc_df = cc_df[cc_df["native"].isin(selected_native_ids)]

                        # Plot using the filtered DataFrame
                        if not filtered_cc_df.empty:
                            filtered_fig_cc, filtered_ax_cc = vis.plot_go_dotplot(
                                filtered_cc_df, category_name="Cellular Component"
                            )
                            filtered_fig_cc.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                            st.pyplot(filtered_fig_cc)
                        else:
                            st.warning("No data available for the selected terms.")

                    # Display the selection table below the plot
                    st.subheader("Select Custom GO Terms")
                    cc_subset_for_selection = cc_df[
                        ["native", "name", "p_value", "q_value", "precision", "recall", "intersections"]
                    ]

                    # Display the dataframe with selections
                    cc_selected_terms = dataframe_with_selections(cc_subset_for_selection, "go_cc_subset_select")

                    # Apply and Reset buttons
                    apply_col, reset_col = st.columns([1, 1])
                    with apply_col:
                        if st.button("Apply Selection", key="cc_apply_selection"):
                            if not cc_selected_terms.empty:
                                st.session_state.cc_selected_terms = cc_selected_terms
                                st.session_state.cc_selection_state = cc_selected_terms.index.tolist()  # Track selected rows
                                st.toast("Selection applied. The plot has been updated!", icon="✅")

                    with reset_col:
                        if st.button("Reset Selection", key="cc_reset_selection"):
                            # Clear selected terms and reset selection state
                            st.session_state.cc_selected_terms = pd.DataFrame()
                            st.session_state.cc_selection_state = []  # Clear the selection state
                            st.toast("Selection reset. Showing default plot.", icon="🔄")



                    
    with go_bp:
        if not enrichment_df.empty and source_dict is not None:
            if "GO:BP" in list(source_dict.keys()):
                bp_df = source_dict["GO:BP"]

                if not bp_df.empty:
                    # Initialize session state for selected terms if not already set
                    if "bp_selected_terms" not in st.session_state:
                        st.session_state.bp_selected_terms = pd.DataFrame()
                    if "bp_selection_state" not in st.session_state:
                        st.session_state.bp_selection_state = []

                    # Render the plot
                    if st.session_state.bp_selected_terms.empty:
                        # Default plot with top 10 GO terms
                        fig_bp, ax_bp = vis.plot_go_dotplot(
                            bp_df.sort_values(by="q_value").head(10),
                            category_name="Biological Process"
                        )
                        fig_bp.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                        st.pyplot(fig_bp)
                    else:
                        # Filter `bp_df` to include only the selected terms
                        selected_native_ids = st.session_state.bp_selected_terms["native"].tolist()
                        filtered_bp_df = bp_df[bp_df["native"].isin(selected_native_ids)]

                        # Plot using the filtered DataFrame
                        if not filtered_bp_df.empty:
                            filtered_fig_bp, filtered_ax_bp = vis.plot_go_dotplot(
                                filtered_bp_df, category_name="Biological Process"
                            )
                            filtered_fig_bp.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                            st.pyplot(filtered_fig_bp)
                        else:
                            st.warning("No data available for the selected terms.")

                    # Display the selection table below the plot
                    st.subheader("Select Custom GO Terms for Biological Process")
                    bp_subset_for_selection = bp_df[
                        ["native", "name", "p_value", "q_value", "precision", "recall", "intersections"]
                    ]
                    bp_selected_terms = dataframe_with_selections(bp_subset_for_selection, "go_bp_subset_select")

                    # Apply and Reset buttons
                    apply_col, reset_col = st.columns([1, 1])
                    with apply_col:
                        if st.button("Apply Selection", key="bp_apply_selection"):
                            if not bp_selected_terms.empty:
                                st.session_state.bp_selected_terms = bp_selected_terms
                                st.session_state.bp_selection_state = bp_selected_terms.index.tolist()
                                st.toast("Selection applied. The plot has been updated!", icon="✅")

                    with reset_col:
                        if st.button("Reset Selection", key="bp_reset_selection"):
                            # Clear selected terms and reset selection state
                            st.session_state.bp_selected_terms = pd.DataFrame()
                            st.session_state.bp_selection_state = []
                            st.toast("Selection reset. Showing default plot.", icon="🔄")
            else:
                st.warning("No Biological Processes were found enriched in GO for this comparison.")

    with go_mf:
        if not enrichment_df.empty and source_dict is not None:
            if "GO:MF" in list(source_dict.keys()):
                mf_df = source_dict["GO:MF"]

                if not mf_df.empty:
                    # Initialize session state for selected terms if not already set
                    if "mf_selected_terms" not in st.session_state:
                        st.session_state.mf_selected_terms = pd.DataFrame()
                    if "mf_selection_state" not in st.session_state:
                        st.session_state.mf_selection_state = []

                    # Render the plot
                    if st.session_state.mf_selected_terms.empty:
                        # Default plot with top 10 GO terms
                        fig_mf, ax_mf = vis.plot_go_dotplot(
                            mf_df.sort_values(by="q_value").head(10),
                            category_name="Molecular Function"
                        )
                        fig_mf.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                        st.pyplot(fig_mf)
                    else:
                        # Filter `mf_df` to include only the selected terms
                        selected_native_ids = st.session_state.mf_selected_terms["native"].tolist()
                        filtered_mf_df = mf_df[mf_df["native"].isin(selected_native_ids)]

                        # Plot using the filtered DataFrame
                        if not filtered_mf_df.empty:
                            filtered_fig_mf, filtered_ax_mf = vis.plot_go_dotplot(
                                filtered_mf_df, category_name="Molecular Function"
                            )
                            filtered_fig_mf.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                            st.pyplot(filtered_fig_mf)
                        else:
                            st.warning("No data available for the selected terms.")

                    # Display the selection table below the plot
                    st.subheader("Select Custom GO Terms for Molecular Function")
                    mf_subset_for_selection = mf_df[
                        ["native", "name", "p_value", "q_value", "precision", "recall", "intersections"]
                    ]
                    mf_selected_terms = dataframe_with_selections(mf_subset_for_selection, "go_mf_subset_select")

                    # Apply and Reset buttons
                    apply_col, reset_col = st.columns([1, 1])
                    with apply_col:
                        if st.button("Apply Selection", key="mf_apply_selection"):
                            if not mf_selected_terms.empty:
                                st.session_state.mf_selected_terms = mf_selected_terms
                                st.session_state.mf_selection_state = mf_selected_terms.index.tolist()
                                st.toast("Selection applied. The plot has been updated!", icon="✅")

                    with reset_col:
                        if st.button("Reset Selection", key="mf_reset_selection"):
                            # Clear selected terms and reset selection state
                            st.session_state.mf_selected_terms = pd.DataFrame()
                            st.session_state.mf_selection_state = []
                            st.toast("Selection reset. Showing default plot.", icon="🔄")
            else:
                st.warning("No Molecular Functions were found enriched for this comparison.")

    with kegg:
        if not enrichment_df.empty and source_dict is not None:
            if "KEGG" in list(source_dict.keys()):
                kegg_df = source_dict["KEGG"]

                if not kegg_df.empty:
                    # Initialize session state for selected terms if not already set
                    if "kegg_selected_terms" not in st.session_state:
                        st.session_state.kegg_selected_terms = pd.DataFrame()
                    if "kegg_selection_state" not in st.session_state:
                        st.session_state.kegg_selection_state = []

                    # Render the plot
                    if st.session_state.kegg_selected_terms.empty:
                        # Default plot with top 10 KEGG pathways
                        fig_kegg, ax_kegg = vis.plot_go_dotplot(
                            kegg_df.sort_values(by="q_value").head(10),
                            category_name="KEGG Pathways"
                        )
                        fig_kegg.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                        st.pyplot(fig_kegg)
                    else:
                        # Filter `kegg_df` to include only the selected terms
                        selected_native_ids = st.session_state.kegg_selected_terms["native"].tolist()
                        filtered_kegg_df = kegg_df[kegg_df["native"].isin(selected_native_ids)]

                        # Plot using the filtered DataFrame
                        if not filtered_kegg_df.empty:
                            filtered_fig_kegg, filtered_ax_kegg = vis.plot_go_dotplot(
                                filtered_kegg_df, category_name="KEGG Pathways"
                            )
                            filtered_fig_kegg.set_size_inches(5, 3)  # Set fixed size: 500px width, 300px height
                            st.pyplot(filtered_fig_kegg)
                        else:
                            st.warning("No data available for the selected terms.")

                    # Display the selection table below the plot
                    st.subheader("Select Custom KEGG Pathways")
                    kegg_subset_for_selection = kegg_df[
                        ["native", "name", "p_value", "q_value", "precision", "recall", "intersections"]
                    ]
                    kegg_selected_terms = dataframe_with_selections(kegg_subset_for_selection, "kegg_subset_select")

                    # Apply and Reset buttons
                    apply_col, reset_col = st.columns([1, 1])
                    with apply_col:
                        if st.button("Apply Selection", key="kegg_apply_selection"):
                            if not kegg_selected_terms.empty:
                                st.session_state.kegg_selected_terms = kegg_selected_terms
                                st.session_state.kegg_selection_state = kegg_selected_terms.index.tolist()
                                st.toast("Selection applied. The plot has been updated!", icon="✅")

                    with reset_col:
                        if st.button("Reset Selection", key="kegg_reset_selection"):
                            # Clear selected terms and reset selection state
                            st.session_state.kegg_selected_terms = pd.DataFrame()
                            st.session_state.kegg_selection_state = []
                            st.toast("Selection reset. Showing default plot.", icon="🔄")
            else:
                st.warning("No KEGG pathways were enriched for this comparison.")
        else:
            st.warning("Enrichment analysis data is not available or no pathways were identified")
    
    with Network_Visualization:
        st.subheader("Protein-GO-KEGG Network Visualization")
    
        # Protein selection
        selected_proteins = st.multiselect(
            "Select Proteins to Display",
            options=dep_list_df["Protein"].tolist(),
            default=dep_list_df["Protein"].head(10).tolist(),
            help="Choose specific proteins to include in the network."
        )
        
        # GO term filtering
        go_categories = ["GO:CC", "GO:BP", "GO:MF"]
        selected_go_categories = st.multiselect(
            "Select GO Categories",
            options=go_categories,
            default=go_categories,
            help="Choose which GO categories to include in the network."
        )
        
        # KEGG pathway filtering
        display_kegg = st.checkbox("Include KEGG Pathways", value=True)

        # Threshold filtering
        p_value_threshold = st.slider(
            "P-value Threshold for Terms",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
            help="Filter GO terms and KEGG pathways by significance."
        )

        # Generate network based on filters
        filtered_source_dict = {
            key: val[val["p_value"] <= p_value_threshold]
            for key, val in source_dict.items()
            if key in selected_go_categories and not val.empty
        }
        filtered_kegg = kegg_df[kegg_df["p_value"] <= p_value_threshold] if display_kegg else pd.DataFrame()

        # Generate the network graph
        network_graph = vis.generate_network_graph(
            dep_list=dep_list_df[dep_list_df["Protein"].isin(selected_proteins)],
            go_terms=filtered_source_dict,
            kegg_terms=filtered_kegg
        )
        
        # Render the graph
        st.components.v1.html(network_graph, height=600)
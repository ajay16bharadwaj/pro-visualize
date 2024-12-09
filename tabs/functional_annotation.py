from utils.decorators import validate_inputs, safe_tab_execution
from config import DEFAULT_COL_DEPLABEL, DEFAULT_FDR_THRESHOLD, DEFAULT_LOG2FC_THRESHOLD, DEFAULT_GENE_NAME_COLUMN 
import streamlit as st
import pandas as pd
from io import BytesIO
from utils.helpers import dataframe_with_selections
from utils.streamlit_caching import cached_get_all_enrichment

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
    ea, go_cc, go_mf, go_bp, kegg = st.tabs(['Comprehensive Enrichment Analysis','GO Cellular Component Encrichment', 'GO Molecular Function', 'GO Biological Process', 'KEGG Biological Pathways'])

    with ea:
        enrichment_df, source_dict =  cached_get_all_enrichment(vis, list(dep_list_df[DEFAULT_GENE_NAME_COLUMN]), organism_input)
        ea_manhattan_plot = vis.plot_manhattan(enrichment_df, category_name="Comprehensive Enrichment Analysis")
        st.plotly_chart(ea_manhattan_plot)

    with go_cc:
        # Column split for table and plot display with a 70-30 ratio.
        go_cc_col1, go_cc_col2 = st.columns([0.7, 0.3])
        if not enrichment_df.empty and source_dict is not None:
            #cc_df = vis.get_go_enrichment(list(dep_list_df['Protein']), go_category="GO:CC", organism=organism_input)
            if "GO:CC" in list(source_dict.keys()):
                cc_df = source_dict['GO:CC']
                
                if not cc_df.empty:
                # Use a checkbox to allow users to select custom GO terms.
                    with go_cc_col2:
                        go_cc_custom_term_select = st.checkbox('Choose custom GO Terms', key='go_cc_custom_term_select')
                        # Display a table for selecting custom GO terms if the checkbox is checked.
                        if go_cc_custom_term_select:
                            cc_selected_terms = dataframe_with_selections(cc_df, "go_cc_df_select")
                            with st.expander("Your selection"):
                                st.write(cc_selected_terms)
                        else:
                            cc_selected_terms = None

                    # Display the plot in the left column.
                    with go_cc_col1:
                        if go_cc_custom_term_select and cc_selected_terms is not None and not cc_selected_terms.empty:
                            # If custom terms are selected, plot based on the selected terms.
                            filtered_fig_cc, filtered_ax_cc = vis.plot_go_dotplot(cc_selected_terms, category_name="Cellular Component")
                            st.pyplot(filtered_fig_cc)
                        else:
                            # Default plot with the top 10 GO terms.
                            fig_cc, ax_cc = vis.plot_go_dotplot(cc_df.sort_values(by='q_value').head(10), category_name="Cellular Component")
                            st.pyplot(fig_cc)

                    #protein information for go term
                    st.write("Select the terms to view the Genes associated")
                    cc_subset_for_selection = cc_df[['native', 'name', 'p_value', 'q_value', 'precision', 'recall','intersections']]
                    cc_selected_terms_for_proteins = dataframe_with_selections(cc_subset_for_selection, "go_cc_protein_select")
                    # Display the list of proteins associated with this term
                    # if not cc_selected_terms_for_proteins.empty:
                    #     st.subheader("Proteins Associated with Selected Term")
                    #     protein_list = cc_selected_terms_for_proteins['intersections'].values[0]  # Access the list of proteins
                    #     protein_df = pd.DataFrame(protein_list, columns=['Protein IDs'])
                    #     st.dataframe(protein_df)

                    if not cc_selected_terms_for_proteins.empty:
                        st.subheader("Genes Associated with Selected Terms")
                        
                        # Prepare data for column-by-column display
                        proteins_dict = {}
                        for _, row in cc_selected_terms_for_proteins.iterrows():
                            term_name = row['name']  # Term name
                            protein_list = row['intersections']  # List of proteins
                            
                            # Add to dictionary, term name as key, and protein list as value
                            proteins_dict[term_name] = protein_list
                        
                        # Create a DataFrame with proteins column by column
                        max_length = max(len(proteins) for proteins in proteins_dict.values())  # Get the maximum protein list length
                        proteins_df = pd.DataFrame({term: pd.Series(proteins) for term, proteins in proteins_dict.items()}, index=range(max_length))
                        
                        # Display the resulting DataFrame
                        st.dataframe(proteins_df)
            else:
                st.warning('No Cellular Componenets were found enriched for this comparison')
                    
        with go_mf:
            go_mf_col1, go_mf_col2 = st.columns([0.7, 0.3])
            if not enrichment_df.empty and source_dict is not None:
                #mf_df = vis.get_go_enrichment(list(dep_list_df['Protein']), go_category="GO:MF", organism=organism_input)
                if "GO:MF" in list(source_dict.keys()):
                    mf_df = source_dict['GO:MF']
                    
                    if not mf_df.empty:
                    # Use a checkbox to allow users to select custom GO terms.
                        with go_mf_col2:
                            go_mf_custom_term_select = st.checkbox('Choose custom GO Terms', key='go_mf_custom_term_select')
                            # Display a table for selecting custom GO terms if the checkbox is checked.
                            if go_mf_custom_term_select:
                                mf_selected_terms = dataframe_with_selections(mf_df, "go_mf_df_select")
                                with st.expander("Your selection"):
                                    st.write(mf_selected_terms)
                            else:
                                mf_selected_terms = None

                        # Display the plot in the left column.
                        with go_mf_col1:
                            if go_mf_custom_term_select and mf_selected_terms is not None and not mf_selected_terms.empty:
                                # If custom terms are selected, plot based on the selected terms.
                                filtered_fig_mf, filtered_ax_mf = vis.plot_go_dotplot(mf_selected_terms, category_name="Cellular Component")
                                st.pyplot(filtered_fig_mf)
                            else:
                                # Default plot with the top 10 GO terms.
                                fig_mf, ax_mf = vis.plot_go_dotplot(mf_df.sort_values(by='q_value').head(10), category_name="Cellular Component")
                                st.pyplot(fig_mf)
                        
                        st.write("Select the terms to view the Genes associated")
                        mf_subset_for_selection = mf_df[['native', 'name', 'p_value', 'q_value', 'precision', 'recall','intersections']]
                        mf_selected_terms_for_proteins = dataframe_with_selections(mf_subset_for_selection, "go_mf_protein_select")
                        # Display the list of proteins associated with this term
                        # if not cc_selected_terms_for_proteins.empty:
                        #     st.subheader("Proteins Associated with Selected Term")
                        #     protein_list = cc_selected_terms_for_proteins['intersections'].values[0]  # Access the list of proteins
                        #     protein_df = pd.DataFrame(protein_list, columns=['Protein IDs'])
                        #     st.dataframe(protein_df)

                        if not mf_selected_terms_for_proteins.empty:
                            st.subheader("Genes Associated with Selected Terms")
                            
                            # Prepare data for column-by-column display
                            proteins_dict = {}
                            for _, row in mf_selected_terms_for_proteins.iterrows():
                                term_name = row['name']  # Term name
                                protein_list = row['intersections']  # List of proteins
                                
                                # Add to dictionary, term name as key, and protein list as value
                                proteins_dict[term_name] = protein_list
                            
                            # Create a DataFrame with proteins column by column
                            max_length = max(len(proteins) for proteins in proteins_dict.values())  # Get the maximum protein list length
                            proteins_df = pd.DataFrame({term: pd.Series(proteins) for term, proteins in proteins_dict.items()}, index=range(max_length))
                            
                            # Display the resulting DataFrame
                            st.dataframe(proteins_df)
                else:
                    st.warning("No Molecular Functions were found enriched with these set of Differentially Expressed Proteins")

        with go_bp:
            go_bp_col1, go_bp_col2 = st.columns([0.7, 0.3])
            if not enrichment_df.empty and source_dict is not None:
                #bp_df = vis.get_go_enrichment(list(dep_list_df['Protein']), go_category="GO:BP", organism=organism_input)
                if "GO:BP" in list(source_dict.keys()):
                    bp_df = source_dict['GO:BP']
                    
                    if not bp_df.empty:
                    # Use a checkbox to allow users to select custom GO terms.
                        with go_bp_col2:
                            go_bp_custom_term_select = st.checkbox('Choose custom GO Terms', key='go_bp_custom_term_select')
                            # Display a table for selecting custom GO terms if the checkbox is checked.
                            if go_bp_custom_term_select:
                                bp_selected_terms = dataframe_with_selections(bp_df, "go_bp_df_select")
                                with st.expander("Your selection"):
                                    st.write(bp_selected_terms)
                            else:
                                bp_selected_terms = None

                        # Display the plot in the left column.
                        with go_bp_col1:
                            if go_bp_custom_term_select and bp_selected_terms is not None and not bp_selected_terms.empty:
                                # If custom terms are selected, plot based on the selected terms.
                                filtered_fig_bp, filtered_ax_bp = vis.plot_go_dotplot(bp_selected_terms, category_name="Cellular Component")
                                st.pyplot(filtered_fig_bp)
                            else:
                                # Default plot with the top 10 GO terms.
                                fig_bp, ax_bp = vis.plot_go_dotplot(bp_df.sort_values(by='q_value').head(10), category_name="Cellular Component")
                                st.pyplot(fig_bp)

                        st.write("Select the terms to view the Genes associated")
                        bp_subset_for_selection = bp_df[['native', 'name', 'p_value', 'q_value', 'precision', 'recall','intersections']]
                        bp_selected_terms_for_proteins = dataframe_with_selections(bp_subset_for_selection, "go_bp_protein_select")
                        # Display the list of proteins associated with this term
                        # if not cc_selected_terms_for_proteins.empty:
                        #     st.subheader("Proteins Associated with Selected Term")
                        #     protein_list = cc_selected_terms_for_proteins['intersections'].values[0]  # Access the list of proteins
                        #     protein_df = pd.DataFrame(protein_list, columns=['Protein IDs'])
                        #     st.dataframe(protein_df)

                        if not bp_selected_terms_for_proteins.empty:
                            st.subheader("Genes Associated with Selected Terms")
                            
                            # Prepare data for column-by-column display
                            proteins_dict = {}
                            for _, row in bp_selected_terms_for_proteins.iterrows():
                                term_name = row['name']  # Term name
                                protein_list = row['intersections']  # List of proteins
                                
                                # Add to dictionary, term name as key, and protein list as value
                                proteins_dict[term_name] = protein_list
                            
                            # Create a DataFrame with proteins column by column
                            max_length = max(len(proteins) for proteins in proteins_dict.values())  # Get the maximum protein list length
                            proteins_df = pd.DataFrame({term: pd.Series(proteins) for term, proteins in proteins_dict.items()}, index=range(max_length))
                            
                            # Display the resulting DataFrame
                            st.dataframe(proteins_df)
                else: 
                    st.warning("No Biological Processes were found enriched in GO for these set of Proteins selected")

            
        with kegg:
            kegg_col1, kegg_col2 = st.columns([0.7, 0.3])
            if not enrichment_df.empty and source_dict is not None:
                #kegg_df = vis.get_go_enrichment(list(dep_list_df['Protein']), go_category="GO:kegg", organism=organism_input)
                #print("source_dict: ", source_dict.keys())
                if "KEGG" in list(source_dict.keys()):
                    kegg_df = source_dict['KEGG']

                    if not kegg_df.empty:
                    
                        # Use a checkbox to allow users to select custom GO terms.
                        with kegg_col2:
                            kegg_custom_term_select = st.checkbox('Choose custom GO Terms', key='kegg_custom_term_select')
                            # Display a table for selecting custom GO terms if the checkbox is checked.
                            if kegg_custom_term_select:
                                kegg_selected_terms = dataframe_with_selections(kegg_df, "kegg_df_select")
                                with st.expander("Your selection"):
                                    st.write(kegg_selected_terms)
                            else:
                                kegg_selected_terms = None

                        # Display the plot in the left column.
                        with kegg_col1:
                            if kegg_custom_term_select and kegg_selected_terms is not None and not kegg_selected_terms.empty:
                                # If custom terms are selected, plot based on the selected terms.
                                filtered_fig_kegg, filtered_ax_kegg = vis.plot_go_dotplot(kegg_selected_terms, category_name="Cellular Component")
                                st.pyplot(filtered_fig_kegg)
                            else:
                                # Default plot with the top 10 GO terms.
                                fig_kegg, ax_kegg = vis.plot_go_dotplot(kegg_df.sort_values(by='q_value').head(10), category_name="Cellular Component")
                                st.pyplot(fig_kegg)

                        st.write("Select the terms to view the Genes associated")
                        kegg_subset_for_selection = kegg_df[['native', 'name', 'p_value', 'q_value', 'precision', 'recall','intersections']]
                        kegg_selected_terms_for_proteins = dataframe_with_selections(kegg_subset_for_selection, "go_kegg_protein_select")
                        # Display the list of proteins associated with this term
                        # if not cc_selected_terms_for_proteins.empty:
                        #     st.subheader("Proteins Associated with Selected Term")
                        #     protein_list = cc_selected_terms_for_proteins['intersections'].values[0]  # Access the list of proteins
                        #     protein_df = pd.DataFrame(protein_list, columns=['Protein IDs'])
                        #     st.dataframe(protein_df)

                        if not kegg_selected_terms_for_proteins.empty:
                            st.subheader("Genes Associated with Selected Terms")
                            
                            # Prepare data for column-by-column display
                            proteins_dict = {}
                            for _, row in kegg_selected_terms_for_proteins.iterrows():
                                term_name = row['name']  # Term name
                                protein_list = row['intersections']  # List of proteins
                                
                                # Add to dictionary, term name as key, and protein list as value
                                proteins_dict[term_name] = protein_list
                            
                            # Create a DataFrame with proteins column by column
                            max_length = max(len(proteins) for proteins in proteins_dict.values())  # Get the maximum protein list length
                            proteins_df = pd.DataFrame({term: pd.Series(proteins) for term, proteins in proteins_dict.items()}, index=range(max_length))
                            
                            # Display the resulting DataFrame
                            st.dataframe(proteins_df)
                else:
                    st.warning("No KEGG pathways were enriched for this set of Differentially Expressed proteins")
            else:
                st.warning("Enrichment analysis data is not available or no pathways were identified")
    
from multiprocessing.managers import ValueProxy
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import scipy.cluster.hierarchy as sch
from scipy.spatial import ConvexHull
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import re
from compute_functions import *
import dash_bio
import requests
from requests.adapters import HTTPAdapter, Retry

class ProteinVisualization:

    COL_DEPLABEL = 'Label2'
    COL_DEPSIGNIF = 'Imputed.FDR'
    color_map = {
            'Down-regulated': '#0072B2',
            'Significant, no change': '#F0E442',
            'Up-regulated': '#D55E00',
            'Non-significant': '#999999'
        }
    
    def __init__(self):
        self.protein_data = None
        self.annotation_info = None
        self.dep_info = None
        self.protein_data_annotated = None
        self.uniprot_info = None
        self.protein_data_for_pca = None
        self.analysis_with_annotation = None
        self.condition_groups = None
        self.protein_data_for_heatmap = None
        self.protein_data_for_violin_plot = None
        self.n_clusters = None
        
    
    def load_protein_data(self, path: str):
        """Load the protein level data from a file"""
        df = pd.read_csv(path, sep='\t')
        self.protein_data = df
        #df.replace([np.inf, -np.inf], np.nan, inplace=True)
        #df.fillna(0, inplace=True)
        #self.protein_data = df.set_index('Protein')

    def load_dep_data(self, path: str):
        """Load the analysis output, DEP data from a file"""
        df = pd.read_csv(path, sep='\t')
        self.dep_info = df

    def load_annotation(self, path: str, SAMPLENAME, GROUPING):
        """Load the annotation information from a file"""
        annotation_df = pd.read_csv(path, sep='\t')
        self.n_clusters = annotation_df[GROUPING].nunique()
        self.annotation_info, self.LIMIT_SAMPLE_COMPARISONS = self.process_annotation_info(annotation_df, SAMPLENAME, GROUPING)
        

    def process_annotation_info(self, annotation_info, SAMPLENAME, GROUPING, LIMIT_SAMPLE_COMPARISONS=False):
        """Process the annotation info and return the processed dataframe."""
        annotation_info = annotation_info.copy()
        annotation_info = annotation_info[[SAMPLENAME, GROUPING]].rename(columns={SAMPLENAME: 'SampleName', GROUPING: 'Group'})
        #self.n_clusters = self.annotation_info['Group'].nunique()
        annotation_info['Group'] = annotation_info['Group'].replace(['', np.nan], 'NOT_ANNOTATED')

        if len(annotation_info['Group'].unique()) > 5 and not LIMIT_SAMPLE_COMPARISONS:
            LIMIT_SAMPLE_COMPARISONS = True
            print("Warning: There are greater than 5 sample groups for this job. --limit_comparisons is turned ON.")
        
        if annotation_info['SampleName'].duplicated().any():
            dups = annotation_info.loc[annotation_info['SampleName'].duplicated(keep=False), 'SampleName'].unique()
            raise ValueError(f"Duplicated sample names found: {', '.join(dups)}")
        
        if annotation_info[['SampleName', 'Group']].isna().any().any():
            raise ValueError(f"No missing values allowed in the annotation file.")
        
        return annotation_info, LIMIT_SAMPLE_COMPARISONS
                

    def preprocess_for_pca(self):
        """Preprocess protein data specifically for PCA or other plots that require normalization."""
        df = self.protein_data.copy()
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.fillna(0, inplace=True)
        #keep only Protein information for all the samples: 
        sample_columns = list(self.annotation_info.SampleName.unique())
        filtered_df = df[['Protein'] + [col for col in sample_columns if col in df.columns]]
        self.protein_data_for_pca = filtered_df.set_index('Protein')
        return self.protein_data_for_pca 

    def preprocess_for_heatmaps(self):
        """Preprocess protein data specifically for heatmaps"""
        df = self.protein_data.copy()
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.fillna(0, inplace=True)
        self.protein_data_for_heatmap = df 
        return self.protein_data_for_heatmap
    
    def preprocess_for_violin_plot(self):
        """Preprocess protein data specifically for heatmaps"""
        df = self.protein_data.copy()
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        df.fillna(0, inplace=True)
        self.protein_data_for_violin_plot = df 
        return self.protein_data_for_violin_plot

    # PCA plot by annotation
    def plot_pca_by_annotation(self, group_column='Group', pc_x=1, pc_y=2):
        """Plot PCA colored by group annotations"""
        # Ensure that the protein data and annotation info have been loaded
        if self.protein_data is None or self.annotation_info is None:
            raise ValueError("Protein data or annotation info not loaded.")

        if self.protein_data_for_pca is None:
            self.preprocess_for_pca()  # Preprocess if not already done

        df = self.protein_data_for_pca
        
        # Step 1: Standardize data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df.T)
    
        # Step 2: Perform PCA
        pca = PCA(n_components=3)
        pca_components = pca.fit_transform(scaled_data)
        pca_df = pd.DataFrame(data=pca_components, columns=[f'PC{i+1}' for i in range(3)])
    
        # Adding jitter to avoid overlap
        jitter = np.random.normal(0, 0.3, size=pca_df.shape)
        pca_df['PC1'] += jitter[:, 0]  # Add jitter to PC1
        pca_df['PC2'] += jitter[:, 1]  # Add jitter to PC2
    
        # Step 3: Merge annotation information
        pca_df = pd.concat([pca_df, self.annotation_info.reset_index()], axis=1)
    
        # Step 4: Create PCA plot colored by annotations
        fig = px.scatter(
            pca_df,
            x=f'PC{pc_x}',
            y=f'PC{pc_y}',
            color=group_column,
            hover_data=['SampleName'],
            title=f'PCA Plot - Colored by {group_column}',
            labels={group_column: 'Groups'},
            text='SampleName',
            color_discrete_map={'Northstar': 'blue', 'Vital': 'green', 'Individual': 'red'}
        )
    
        # Dynamically adjust the axis scale and the labels positioning
        fig.update_traces(textposition='top center')
        
        fig.update_layout(
            xaxis_title=f'PC{pc_x} ({pca.explained_variance_ratio_[pc_x-1]*100:.1f}% Variance)',
            yaxis_title=f'PC{pc_y} ({pca.explained_variance_ratio_[pc_y-1]*100:.1f}% Variance)',
            width=800,
            height=600,
            margin=dict(l=0, r=0, t=40, b=40),
            showlegend=True,
        )
        
        # Adjust marker size
        fig.update_traces(marker=dict(size=8))
    
        return fig

    # PCA plot with clusters
    def plot_pca_with_clusters_plotly(self):
        """Plot PCA with clusters using Plotly and Convex Hulls"""
        from scipy.spatial import ConvexHull
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        from sklearn.cluster import KMeans
        from matplotlib.patches import Ellipse
        
        # Ensure protein data is loaded
        if self.protein_data is None:
            raise ValueError("Protein data not loaded.")

        if self.protein_data_for_pca is None:
            self.preprocess_for_pca()  # Preprocess if not already done

        df = self.protein_data_for_pca
        
        # Standardize the data
        #df_scaled = StandardScaler().fit_transform(df.values.T)
        #print(df.dtypes)
        #print(df.head())

        # Standardize the data and convert to float64
        df_scaled = StandardScaler().fit_transform(df.values.T.astype(np.float64))
    
        # Perform PCA
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(df_scaled)
    
        # Clustering
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)
        clusters = kmeans.fit_predict(principal_components)
    
        # Create a dataframe for PCA results and cluster labels
        pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
        pca_df['Cluster'] = clusters
        pca_df['Sample'] = df.columns # type: ignore
    
        # Initialize plot
        fig = go.Figure()
    
        colors = ['#FF6347', '#4682B4', '#32CD32']  # Define colors for clusters
    
        for cluster in range(self.n_clusters): # type: ignore
            # Filter points of the current cluster
            cluster_points = pca_df[pca_df['Cluster'] == cluster]
    
            # Plot points for the cluster
            fig.add_trace(go.Scatter(
                x=cluster_points['PC1'], y=cluster_points['PC2'],
                mode='markers+text',
                text=cluster_points['Sample'],
                marker=dict(size=10, color=colors[cluster]),
                name=f'Cluster {cluster}'
            ))
    
            # Convex Hull
            points = cluster_points[['PC1', 'PC2']].values
            if len(points) > 2:  # ConvexHull requires at least 3 points
                hull = ConvexHull(points)
                hull_points = np.append(hull.vertices, hull.vertices[0])  # Close the hull
                fig.add_trace(go.Scatter(
                    x=points[hull_points, 0], y=points[hull_points, 1],
                    fill='toself',
                    fillcolor=colors[cluster],
                    opacity=0.2,
                    line=dict(color=colors[cluster]),
                    showlegend=False
                ))
    
        # Update layout
        fig.update_layout(
            title="PCA - Colored by Clusters (PC1 vs PC2)",
            xaxis_title=f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)',
            yaxis_title=f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)',
            legend_title="Cluster",
            autosize=False,
            width=1200,  # Adjust width to make the figure longer
            height=600,  # Adjust height as per your requirement
        )
    
        return fig

    def plot_vertical_dendrogram(self, method='ward', figsize=(8, 10), label_size=10):
        """Plot vertical dendrogram using hierarchical clustering."""
    
        # Ensure protein data is loaded
        if self.protein_data is None:
            raise ValueError("Protein data not loaded.")

        if self.protein_data_for_pca is None:
            self.preprocess_for_pca()  # Preprocess if not already done

        df = self.protein_data_for_pca
        
        # Standardize the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df.T)
    
        # Generate the linkage matrix
        Z = sch.linkage(scaled_data, method=method)
    
        # Create the figure for the dendrogram
        fig, ax = plt.subplots(figsize=figsize)
        plt.title('Hierarchical Clustering Dendrogram (Vertical)', fontsize=14)
        
        # Plot vertical dendrogram
        sch.dendrogram(Z, labels=df.columns, orientation='right', leaf_font_size=label_size, ax=ax)
    
        # Adjust label sizes and spacing
        plt.xticks(fontsize=label_size)
        plt.yticks(fontsize=label_size)
        
        plt.xlabel("Distance")
        plt.ylabel("Samples")
        plt.tight_layout()

        plt.close(fig)
    
        return fig

    def create_cluster_assignment_table(self):
        """Create a table with cluster assignments for each sample"""
        from sklearn.cluster import AgglomerativeClustering
        
        # Ensure protein data is loaded
        if self.protein_data is None:
            raise ValueError("Protein data not loaded.")

        if self.protein_data_for_pca is None:
            self.preprocess_for_pca()  # Preprocess if not already done

        df = self.protein_data_for_pca
        
        # Standardize data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df.T)
    
        # Perform clustering
        clustering = AgglomerativeClustering(n_clusters=self.n_clusters)
        cluster_labels = clustering.fit_predict(scaled_data)
    
        # Create a DataFrame for the assignments
        cluster_assignment_df = pd.DataFrame({
            'SampleName': df.columns,
            'Cluster': cluster_labels
        })
    
        return cluster_assignment_df

    def get_uniprot_info(self):
        """Fetch UniProt information for the proteins in the dataset."""
        if self.uniprot_info is None:
           
            proteins =   list(set(self.protein_data['Protein'].tolist()))# type: ignore # Get unique proteins from the index
            protein_names = []

            # Remove fasta descriptions if present and get UniProt IDs
            for protein_name in proteins:
                pt = self._get_uniprot_ids(protein_name)
                protein_names.append(pt)

            # Split the protein list into chunks
            chunk_size = 130
            if len(protein_names) > chunk_size:
                split_list = [protein_names[i:i+chunk_size] for i in range(0, len(protein_names), chunk_size)]
            else:
                split_list = [protein_names]

            # Define the columns to retrieve from UniProt
            cols = ["accession", "gene_primary", "protein_name", "cc_tissue_specificity", "go_p", "go_c", "go_f", "cc_subcellular_location"]
            cols = ",".join(cols)

            # Get UniProt entries
            uniprot_df = self.get_uniprot_entries(protein_names, split_list, cols)
            self.uniprot_info = uniprot_df

        

    def _get_uniprot_ids(self, row):
        """Helper function to extract UniProt IDs from protein names"""
        REGEX_UNIPROTID = re.compile("iRT|([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?")
        result = re.search(REGEX_UNIPROTID, row)
        if result is None:
            pt = row
        else:
            pt = result.group(0)
        return pt

    def get_uniprot_entries(self, protein_names, split_list, cols):
        """Fetch entries from UniProt using API."""
        def get_next_link(headers):
            if "Link" in headers:
                match = re_next_link.match(headers["Link"])
                if match:
                    return match.group(1)

        def get_batch(batch_url):
            while batch_url:
                response = session.get(batch_url)
                response.raise_for_status()
                yield response
                batch_url = get_next_link(response.headers)

        # Set up session for requests with retries
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))

        # Store entries
        entries = []

        # Iterate through split protein lists to make API requests
        for protein_list in split_list:
            query_proteins = "+OR+".join(protein_list)
            url = f'https://rest.uniprot.org/uniprotkb/search?fields={cols}&size=300&format=tsv&query={query_proteins}'

            for batch in get_batch(url):
                for line in batch.text.splitlines()[1:]:
                    line = line.split("\t")
                    if line[0] not in protein_names:
                        continue
                    if line not in entries:
                        entries.append(line)

        if len(entries) != 0:
            # Define the column headers for the UniProt DataFrame
            df_cols = ["ProteinIds", "Gene Name", "Protein Description", "Tissue Specificity", 
                       "Gene ontology (biological process)", "Gene ontology (cellular component)", 
                       "Gene ontology (molecular function)", "Subcellular Location[CC]"]
            return pd.DataFrame(entries, columns=df_cols)
        else:
            print("Couldn't fetch UniProt Annotation at the moment.")
            return pd.DataFrame()

    def merge_uniprot2proteome(self, df):
        """Merge DEP data with UniProt information."""
        if self.uniprot_info is None:
            raise ValueError("UniProt information is not loaded.")
        df['ProteinIds'] = df['Protein'].apply(self._get_uniprot_ids)
        merged_df = pd.merge(self.uniprot_info, df, on="ProteinIds", how="right")
        return merged_df
    
    def merge_dep_ptlevel(self, dep_df, pt_level):
        df2_selected = pt_level[['Protein', 'Gene Name', 'Protein Description']]
        merged_df = pd.merge(dep_df, df2_selected, on="Protein", how="left")
        return merged_df

    def volcano_preprocess(self, COL_DEPLABEL):
        """Preprocess for volcano plot by grouping based on DEP label."""
        if self.dep_info is None:
            raise ValueError("DEP information is not loaded or Uniprot Annotation is non existant")
        
        if self.protein_data_annotated is None:
            raise ValueError("Protein Data is not annotated, please upload annotated protein information")
        
        if self.analysis_with_annotation is None:        
            self.analysis_with_annotation = self.merge_dep_ptlevel(self.dep_info, self.protein_data_annotated)
        
        condition_level = sorted(self.analysis_with_annotation[COL_DEPLABEL].unique())
        self.analysis_with_annotation[COL_DEPLABEL] = pd.Categorical(self.analysis_with_annotation[COL_DEPLABEL], categories=condition_level, ordered=True)
        condition_groups = {condition: group for condition, group in self.analysis_with_annotation.groupby(COL_DEPLABEL)}
        return condition_groups
    
    def plot_volcano(self, df, FC, p_val):
        """Generate volcano plot."""
        if self.protein_data_annotated is None: 
            raise ValueError("Annotated Protein data is required to display Gene and Protein Description information")
        
    

        df['sig'] = -np.log10(df[self.COL_DEPSIGNIF])

        df['label'] = 'Non-significant'
        df.loc[df[self.COL_DEPSIGNIF] <= p_val, 'label'] = 'Significant, no change'
        df.loc[(df[self.COL_DEPSIGNIF] <= p_val) & (df['log2FC'] > FC), 'label'] = 'Up-regulated'
        df.loc[(df[self.COL_DEPSIGNIF] <= p_val) & (df['log2FC'] < -FC), 'label'] = 'Down-regulated'
    
        fig = px.scatter(df, x='log2FC', y='sig', color='label', color_discrete_map=self.color_map, hover_data=['Gene Name', 'Protein Description'])
            
            # Customizing the layout
        fig.update_layout(title=f'Volcano Plot', xaxis_title="log2 fold change", yaxis_title="-log10(FDR)", legend_title_text='Category')
            
            # Adding horizontal and vertical lines for cutoffs
        fig.add_shape(type="line", x0=-FC, y0=0, x1=-FC, y1=df['sig'].max(), line=dict(color="RoyalBlue", width=2, dash="dot"))
        fig.add_shape(type="line", x0=FC, y0=0, x1=FC, y1=df['sig'].max(), line=dict(color="RoyalBlue", width=2, dash="dot"))
        fig.add_shape(type="line", x0=df['log2FC'].min(), y0=-np.log10(p_val), x1=df['log2FC'].max(), y1=-np.log10(p_val),line=dict(color="Red", width=2, dash="dot"))
    
        return fig

    def plot_protein_overlap(self):
        """
        Plot the protein overlap showing the number of proteins identified across samples.
        Args:
            protein_data: DataFrame with proteins as rows and samples as columns (with numeric values).
        Returns:
            fig: Plotly figure of the protein overlap across samples.
        """
        protein_data = self.protein_data.copy()
    
        # Remove any non-numeric columns (e.g., protein names or descriptions)
        numeric_data = protein_data.select_dtypes(include=[float, int])
    
        # Count how many times each protein is identified across samples (non-zero values)
        protein_counts = (numeric_data > 0).sum(axis=1)
    
        # Count how many proteins appear in exactly 'X' samples
        overlap_counts = protein_counts.value_counts().sort_index()
    
        # Plot: Protein overlap across all samples (bar plot)
        fig = go.Figure()
    
        fig.add_trace(go.Bar(
            x=overlap_counts.index,
            y=overlap_counts.values,
            marker=dict(color='black'),
            name="Protein Overlap"
        ))
    
        fig.update_layout(
            title="Protein Overlap Across Samples",
            xaxis_title="Identified in X samples",
            yaxis_title="Number of shared proteins",
            width=600,
            height=400
        )
    
        return fig
    
    def plot_proteins_per_sample(self, group_column='Group'): # type: ignore
        """
        Creates a bar plot of the number of identified proteins per sample, dynamically colored by groups.
        """
        # Preprocessing to count proteins per sample (non-NaN values)
        protein_data  = self.protein_data.copy()
        annotation_info = self.annotation_info.copy()
        protein_counts_per_sample = protein_data.iloc[:, 1:].notna().sum()

        # Convert to DataFrame and merge with annotation information
        protein_counts_df = pd.DataFrame({
            'SampleName': protein_data.columns[1:],  # Sample names as columns (exclude the Protein column)
            'ProteinCount': protein_counts_per_sample.values
        })
        
        # Merge with annotation data to include groups
        merged_df = pd.merge(protein_counts_df, annotation_info, how='left', on='SampleName')

        # Dynamically assign colors based on unique groups in the annotation
        unique_groups = merged_df[group_column].unique()
        colors = px.colors.qualitative.Plotly[:len(unique_groups)]  # Select colors dynamically
        color_discrete_map = {group: color for group, color in zip(unique_groups, colors)}

        # Calculate the dynamic threshold for the horizontal line
        max_protein_count = merged_df['ProteinCount'].max()
        dynamic_threshold = max_protein_count + 100  # Add 100 to the highest value

        # Plotting
        fig = px.bar(
            merged_df,
            x='SampleName',
            y='ProteinCount',
            color=group_column,
            title="Identified Proteins per Sample",
            labels={'ProteinCount': 'Number of Proteins', 'SampleName': 'Sample Name'},
            color_discrete_map=color_discrete_map
        )

        # Add a dynamic horizontal line based on the max protein count
        fig.add_shape(
            type="line",
            x0=0,
            x1=1,
            y0=dynamic_threshold,  # Dynamic threshold
            y1=dynamic_threshold,
            line=dict(color="black", width=2, dash="dash"),
            xref='paper',
            yref='y'
        )

        # Update layout for better label placement and readability
        fig.update_layout(
            xaxis_tickangle=-45,  # Rotate x-axis labels for better readability
            xaxis_tickfont=dict(size=12),  # Adjust font size
            yaxis_title="Number of Proteins",
            xaxis_title="Sample Name",
            height=500,
            width=900,
            margin=dict(l=50, r=50, t=100, b=150),  # Adjust bottom margin for larger labels
            legend_title="Group"
        )

        return fig
    
    def plot_violin_with_subplots(self, protein_list):
        from plotly.subplots import make_subplots
        """
        Generate a violin plot for one or more proteins with subplots, dynamically colored by groups.
        Args:
            protein_list: List of proteins to generate the plot for (can be a single protein or a list).
        Returns:
            fig: Plotly figure of the violin plot(s) in subplots.
        """
        # Ensure protein data and annotation info are loaded
        if self.protein_data is None or self.annotation_info is None:
            raise ValueError("Protein data or annotation info not loaded.")
        
        # Ensure the protein_list is a list
        if isinstance(protein_list, str):
            protein_list = [protein_list]
        
        df = self.protein_data.copy()
        # Filter protein data for the specified proteins
        #ensure only sample columns are present. 
        sample_columns = list(self.annotation_info.SampleName.unique())
        filtered_df = df[['Protein'] + [col for col in sample_columns if col in df.columns]]


        #need to set the index and scale it before plotting. 

        df_filtered = filtered_df[filtered_df['Protein'].isin(protein_list)]
        if df_filtered.empty:
            raise ValueError("None of the specified proteins were found in the protein data.")
        
        # Melt the DataFrame (wide to long format)
        df_long = pd.melt(df_filtered, id_vars='Protein', var_name='Sample', value_name='Expression')
        
        # Merge with annotation info to get the group labels for each sample
        df_long = pd.merge(df_long, self.annotation_info, how='left', left_on='Sample', right_on='SampleName')
        
        # Dynamically assign colors based on unique groups in the annotation
        unique_groups = df_long['Group'].unique()
        colors = px.colors.qualitative.Plotly[:len(unique_groups)]  # Select colors dynamically
        color_discrete_map = {group: color for group, color in zip(unique_groups, colors)}
        
        # Set up the subplot structure
        num_proteins = len(protein_list)
        cols = 3  # Define the number of columns for the subplots
        rows = (num_proteins // cols) + 1  # Calculate the number of rows based on the number of proteins

        fig = make_subplots(rows=rows, cols=cols, subplot_titles=protein_list)

        # Track which groups have already been added to the legend
        shown_groups = set()

        # Add violin plots for each protein in subplots
        for i, protein in enumerate(protein_list):
            protein_data = df_long[df_long['Protein'] == protein]
            row = (i // cols) + 1
            col = (i % cols) + 1
            
            # Get min and max values for y-axis dynamically
            y_min = protein_data['Expression'].min()
            y_max = protein_data['Expression'].max()

            for group in protein_data['Group'].unique():
                group_data = protein_data[protein_data['Group'] == group]

                # Only show the legend for the first occurrence of each group
                show_legend = group not in shown_groups
                shown_groups.add(group)

                fig.add_trace(
                    go.Violin(
                        x=group_data['Group'],
                        y=group_data['Expression'],
                        name=group,
                        box_visible=True,  # Show the boxplot within the violin plot
                        meanline_visible=True,  # Show the mean line
                        pointpos=0,  # Position of the scatter points
                        jitter=0.05,  # Jitter for scatter points
                        scalegroup=protein,
                        side='both',  # Plot on both sides of the axis for wider violin plots
                        fillcolor=color_discrete_map[group],  # Apply the group-specific color
                        line_color="black",  # Color for the box outline
                        width=0.8,  # Make the violin plots wider
                        box_fillcolor="white",  # White box for the middle
                        showlegend=show_legend  # Only show legend for the first occurrence
                    ),
                    row=row,
                    col=col
                )

            # Set dynamic y-axis range for each subplot and apply scientific notation (1e6, 2e6, etc.)
            fig.update_yaxes(
                range=[y_min - y_min * 0.1, y_max + y_max * 0.1],
                row=row,
                col=col,
                exponentformat="e",  # Scientific notation
                showexponent='all'
            )

        # Add the legend manually to show color associations
        fig.update_layout(
            title=f"Violin Plot - Top Differentially Expressed Proteins",
            height=400 * rows,  # Adjust the height based on the number of rows
            width=1200,  # Adjust the width of the plot
            showlegend=True,
            legend=dict(
                title="Group",
                itemsizing='constant',
                font=dict(size=12),
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="right",
                x=1.2
            ),
            template='plotly_white',
            font=dict(size=12),  # Adjust font size
            margin=dict(l=50, r=50, t=100, b=50)  # Tighten margins
        )

        return fig

    def plot_heatmap(self, df):

        from sklearn.preprocessing import MinMaxScaler

        #ensure only numeric values from samples get passed on to the plot. 
        # Get the list of sample names from annotation information
        sample_columns = list(self.annotation_info.SampleName.unique())
        filtered_df = df[['Protein'] + [col for col in sample_columns if col in df.columns]]


        #need to set the index and scale it before plotting. 
        filtered_df = filtered_df.set_index('Protein')
        #df = df.drop('Protein', axis=1)

        scaler = MinMaxScaler()
        df1_norm = pd.DataFrame(scaler.fit_transform(filtered_df), columns=filtered_df.columns, index=filtered_df.index)


        fig = dash_bio.Clustergram(
        data=df1_norm,
        column_labels=list(df1_norm.columns.values),
        row_labels=list(df1_norm.index),
        height=600,
        width=600
        )
        
        fig.update_layout(title=f'Heatmap', xaxis_title="samples", yaxis_title="Proteins") # type: ignore

        return fig  
    
    def plot_correlation_matrix(self):
        """Plot a correlation matrix using Spearman correlation from raw protein data and adjust layout for better visualization."""
        # Ensure protein data is loaded
        if self.protein_data is None:
            raise ValueError("Protein data not loaded.")

        # Use the raw protein data, setting 'Protein' as the index
        df = self.protein_data.copy()
        sample_columns = list(self.annotation_info.SampleName.unique())
        filtered_df = df[['Protein'] + [col for col in sample_columns if col in df.columns]]
        filtered_df = filtered_df.set_index('Protein')

        # Replace infinite values and fill NaNs with zeros
        filtered_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        filtered_df.fillna(0, inplace=True)

        # Calculate Spearman correlation
        corr_matrix = filtered_df.corr(method='spearman')

        # Create the heatmap using Plotly
        fig = go.Figure(
            data=go.Heatmap(
                z=corr_matrix.values,
                x=corr_matrix.columns,
                y=corr_matrix.index,
                colorscale='RdBu',
                zmin=-1,
                zmax=1,
                colorbar=dict(title='Correlation', titleside='right', tickvals=[-1, -0.5, 0, 0.5, 1]),
            )
        )

        # Update layout for better readability
        fig.update_layout(
            title='Correlation Matrix - Spearman Correlation',
            xaxis=dict(
                title='Samples',
                tickangle=45,
                tickfont=dict(size=8),
                automargin=True,
                side='bottom'
            ),
            yaxis=dict(
                title='Samples',
                tickfont=dict(size=8),
                automargin=True,
            ),
            width=800,
            height=800,
            margin=dict(l=100, r=20, t=100, b=100),
            title_font=dict(size=14)
        )
        
        # Set aspect ratio to make it look square
        fig.update_yaxes(scaleanchor="x", scaleratio=1)

        return fig
    
    def plot_umap(self):
        import umap
        """
        Generates an improved UMAP plot with dynamic color assignment based on sample groups.
        
        Parameters:
        - self: Instance of the ProteinVisualization class.
        
        Returns:
        - fig: Plotly figure object for UMAP projection.
        """
        # Ensure that protein data and annotation info are loaded
        if self.protein_data is None or self.annotation_info is None:
            raise ValueError("Protein data or annotation information not loaded.")

        # Preprocess protein data for UMAP
        df = self.preprocess_for_pca()  # You can reuse the PCA preprocessing method for UMAP
        df = df.T  # Transpose to have samples as rows and proteins as columns

        # Ensure that the dimensions match between annotation and data
        df = df.loc[self.annotation_info['SampleName']]  # Ensure samples match the annotation info
        
        # Perform UMAP on the preprocessed protein data
        reducer = umap.UMAP(random_state=42, n_neighbors=15, min_dist=0.1)
        umap_embedding = reducer.fit_transform(df)

        # Add UMAP projection columns to the annotation info DataFrame
        self.annotation_info['UMAP1'] = umap_embedding[:, 0]
        self.annotation_info['UMAP2'] = umap_embedding[:, 1]

        # Get unique groups dynamically from the annotation file and create color mapping
        unique_groups = self.annotation_info['Group'].unique()
        color_palette = px.colors.qualitative.Plotly  # Dynamic color palette from Plotly
        color_discrete_map = {group: color_palette[i % len(color_palette)] for i, group in enumerate(unique_groups)}

        # Create UMAP plot using Plotly Express
        fig = px.scatter(
            self.annotation_info,
            x='UMAP1',
            y='UMAP2',
            color='Group',
            hover_data=['SampleName'],
            title="UMAP Projection",
            labels={'UMAP1': 'UMAP1', 'UMAP2': 'UMAP2'},
            color_discrete_map=color_discrete_map
        )

        # Update layout to improve aesthetics
        fig.update_traces(marker=dict(size=10))
        fig.update_layout(
            height=600,
            width=800,
            legend_title_text="Groups",
            title_x=0.5,
            title_font=dict(size=20),
            font=dict(size=12),
            margin=dict(l=40, r=40, t=60, b=40)
        )

        return fig
        

    def plot_tsne(self, group_column='Group', n_components=2, perplexity=10):
        from sklearn.preprocessing import StandardScaler
        from sklearn.manifold import TSNE
        """Generate a t-SNE plot for visualizing protein data."""
        # Ensure that the protein data and annotation info have been loaded
        if self.protein_data is None or self.annotation_info is None:
            raise ValueError("Protein data or annotation info not loaded.")

        if self.protein_data_for_pca is None:
            self.preprocess_for_pca()  # Reusing preprocessing from PCA

        df = self.protein_data_for_pca
        
        n_samples = df.T.shape[0]  # Number of samples
    
        # Adjust the perplexity if it's higher than n_samples - 1
        if perplexity >= n_samples:
            perplexity = max(5, n_samples // 2)
        
        # Standardize the data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df.T)
        
        # Perform t-SNE
        tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
        tsne_components = tsne.fit_transform(scaled_data)
        tsne_df = pd.DataFrame(data=tsne_components, columns=[f'TSNE{i+1}' for i in range(n_components)])
        
        # Step 3: Merge annotation information
        tsne_df = pd.concat([tsne_df, self.annotation_info.reset_index()], axis=1)
        
        # Step 4: Create t-SNE plot colored by annotations
        fig = px.scatter(
            tsne_df,
            x='TSNE1',
            y='TSNE2',
            color=group_column,
            hover_data=['SampleName'],
            title=f't-SNE Plot - Colored by {group_column}',
            labels={group_column: 'Groups'},
            text='SampleName'
        )
        
        # Dynamically adjust the axis scale and the labels positioning
        fig.update_traces(textposition='top center')
        
        fig.update_layout(
            xaxis_title='t-SNE Component 1',
            yaxis_title='t-SNE Component 2',
            width=800,
            height=600,
            margin=dict(l=0, r=0, t=40, b=40),
            showlegend=True,
        )
        
        # Adjust marker size
        fig.update_traces(marker=dict(size=8))
        
        return fig
    
    def plot_intensity_density(self):
        """
        Create density plots of protein intensities grouped by sample groups.
        Assumes that the annotation information includes a 'Group' column
        and that the intensity data is stored in 'protein_data'.
        """
        # Ensure protein data and annotation info are loaded
        if self.protein_data is None or self.annotation_info is None:
            raise ValueError("Protein data or annotation info not loaded.")
        
        # Merge protein data with annotation info
        df = self.protein_data.copy()
        #ensure only sample columns are present. 
        sample_columns = list(self.annotation_info.SampleName.unique())
        filtered_df = df[['Protein'] + [col for col in sample_columns if col in df.columns]]


        filtered_df = pd.melt(filtered_df, id_vars=['Protein'], var_name='Sample', value_name='Intensity')
        
        # Log transformation to make the data more normally distributed
        filtered_df['log10(Intensity)'] = np.log10(filtered_df['Intensity'].replace(0, np.nan))
        
        # Merge with annotation info to add grouping
        filtered_df = filtered_df.merge(self.annotation_info, left_on='Sample', right_on='SampleName', how='left')
        
        # Create the density plot using plotly express
        fig = px.histogram(
            filtered_df,
            x='log10(Intensity)',
            color='Sample',
            facet_col='Group',
            histnorm='density',
            nbins=100,
            opacity=0.7,
            title='Density Plot of Protein Intensities by Group',
            labels={'log10(Intensity)': 'log10(Intensity)', 'Group': 'Group'},
            height=600,
            width=1000
        )
        
        # Update layout for better visibility
        fig.update_layout(
            title_text='Density Plot of Protein Intensities by Group',
            legend_title='Sample',
            margin=dict(t=50, l=50, r=50, b=50),
            xaxis_title='log10(Intensity)',
            yaxis_title='Density'
        )
        
        # Adjust facet titles for better presentation
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        
        return fig
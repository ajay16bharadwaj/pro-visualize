# **Release Pro-Visualize v1.0.0: Major Updates and New Features**

We are excited to announce the release of Pro-Visualize **v1.0.0**, which introduces a range of new features, enhancements, and bug fixes to enhance the platform's functionality, usability, and visualizations. This release represents a significant step forward in delivering advanced tools for proteomics data analysis and visualization.

---

## **Highlights of this Release**

### **1. Visualization Improvements**
- **Protein Rank Order Diagram**: A new visualization feature for ranking proteins based on quantitative data.
- **Venn Diagram Visualizations**:
  - Visualize overlapping protein groups.
  - Compare proteins across different datasets.
- **Volcano Plot Enhancements**:
  - Added cutoff information to improve clarity and usability.
- **PCA Visualizations**:
  - Dynamic coloring support for better clustering analysis.
  - Simplified visuals by removing redundant text labels.

### **2. Functional Annotation and Pathway Analysis**
- Introduced **functional annotation** features, now operating at the gene level.
- Integrated **KEGG** and **GO_BP** pathway analysis:
  - Beta functionality to select and display proteins associated with specific terms.

### **3. Customization and Interactivity**
- Enhanced **heatmap** and **violin plots** with options for custom protein subset selection.
- Added functionality to display selected proteins in the **KEGG** and **GO_BP** tabs.

---

## **Bug Fixes and Refinements**
- Fixed **Uniprot annotation** issues; annotations now operate at the protein level.
- Resolved indexing issues in **comparison groups**.
- Improved error handling for smoother operation.
- Removed redundant print statements for cleaner logs.

---

## **Miscellaneous Updates**
- Updated `requirements.txt` to include new dependencies.
- Enhanced `.gitignore` to improve repository hygiene.

---

## **Getting Started**
To use this release:
1. Pull the latest changes from the `main` branch or download the release package.
2. Ensure all dependencies are installed by running:
   ```bash
   pip install -r requirements.txt


streamlit run streamlit_app.py

import streamlit as st
import pandas as pd
from config import DEFAULT_VOLCANO_CONFIG

## should I move this to compute or helper functions? 
def dataframe_with_selections(df, custom_key_name):
    """Provide an interactive selection option in a Streamlit data editor."""
    df_with_selections = df.copy()
    df_with_selections.insert(0, "Select", False)

    # Get dataframe row-selections from user with st.data_editor 
    # Interactive selection
    edited_df = st.data_editor(
        df_with_selections,
        hide_index=True,
        column_config={"Select": st.column_config.CheckboxColumn(required=True)},
        disabled=df.columns,
        key=custom_key_name
    )

    # Filter the dataframe using the temporary column, then drop the column
    # Return selected rows
    selected_rows = edited_df[edited_df.Select]
    return selected_rows.drop('Select', axis=1)

def get_user_volcano_config():
    user_config = DEFAULT_VOLCANO_CONFIG.copy()
    
    if st.checkbox("Customize Volcano Plot"):
        user_config["title"] = st.text_input("Plot Title", user_config["title"])
        user_config["x_label"] = st.text_input("X-axis Label", user_config["x_label"])
        user_config["y_label"] = st.text_input("Y-axis Label", user_config["y_label"])
        user_config["label_size"] = st.slider("Label Font Size", 10, 30, user_config["label_size"])
        
        # Color customization
        user_config["colors"]["Up-regulated"] = st.color_picker("Up-regulated Color", user_config["colors"]["Up-regulated"])
        user_config["colors"]["Down-regulated"] = st.color_picker("Down-regulated Color", user_config["colors"]["Down-regulated"])
        user_config["colors"]["Non-significant"] = st.color_picker("Non-significant Color", user_config["colors"]["Non-significant"])
        user_config["colors"]["Significant, no change"] = st.color_picker("Significant, no change Color", user_config["colors"]["Significant, no change"])
        
        # Threshold line toggle
        user_config["threshold_lines"] = st.checkbox("Show Threshold Lines", user_config["threshold_lines"])

    return user_config

import streamlit as st
import pandas as pd


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



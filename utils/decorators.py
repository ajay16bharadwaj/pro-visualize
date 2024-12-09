import streamlit as st
from functools import wraps

def validate_inputs(analysis_status=False, protein_status=False, annotation_status=False):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Collect missing inputs
            missing_inputs = []

            if analysis_status and not kwargs.get('analysis_status', False):
                missing_inputs.append("Analysis input")
            if annotation_status and not kwargs.get('annotation_status', False):
                missing_inputs.append("Annotation input")
            if protein_status and not kwargs.get('protein_status', False):
                missing_inputs.append("Protein level input")

            # If any inputs are missing, show a combined warning
            if missing_inputs:
                st.warning(f"{', '.join(missing_inputs)} not provided. Please upload the required files.")
                return

            # Call the original function if all checks pass
            return func(*args, **kwargs)
        return wrapper
    return decorator






##decorator function for execption handling
def safe_tab_execution(tab_name):
    """Decorator to handle errors within a tab."""
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                func(*args, **kwargs)
            except Exception as e:
                st.error(f"An error occurred in the '{tab_name}' tab: {str(e)}")
        return wrapper
    return decorator
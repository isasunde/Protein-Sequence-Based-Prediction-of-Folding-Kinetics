import numpy as np
import pandas as pd

def data_filter(
        data: pd.DataFrame, 
        columns_of_interest: list
        ):
    """
    Filter the dataset to include only the specified columns and remove rows with null values.

    Parameters:
    data (pd.DataFrame): The input dataset.
    columns_of_interest (list): List of column names to retain.

    Returns:
    pd.DataFrame: The filtered dataset.
    """
    # Select relevant columns
    data = data[columns_of_interest]
    
    # Remove rows with null values in any of the specified columns
    data = data.dropna(subset=columns_of_interest)
    
    return data
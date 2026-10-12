import json
import requests
import numpy as np
import pandas as pd

from datetime import datetime, timezone

def get_endpoint(
    id: str,
    API_URL: str,
    endpoint: str = 'uniprot',
    params: dict | None = None,
) -> tuple[str, (dict | list)]:
    """Send API request in the K-pro format"""
    # Format API request
    url = f"{API_URL}/{endpoint}/{id}"

    # Send API request
    try:
        response = requests.get(
            url,
            params=params,
            timeout=30,
        )

        response.raise_for_status()

    # Handle potential API errors
    except requests.exceptions.HTTPError as http_error_message:
        print (f"[HTTP ERROR]: {http_error_message}")
    
    except requests.exceptions.ConnectionError as connection_error_message:
        print (f"[CONNECTION ERROR]: {connection_error_message}")
    
    except requests.exceptions.Timeout as timeout_error_message:
        print (f"[TIMEOUT ERROR]: {timeout_error_message}")
    
    except requests.exceptions.RequestException as other_error_message:
        print (f"[UNKNOWN ERROR]: {other_error_message}")

    return url, response

def save_raw_response(
    data: list,
    url: str,
    endpoint: str,
    filename: str,
    RAW_DIR: str,
):
    """Save an API response together with request metadata."""
    # Check if the given directory exists or create it
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    # Collect metadata and API response
    payload = {
        "metadata": {
            "extracted_at": datetime.now(timezone.utc).isoformat(),
            "url": url,
            "endpoint": endpoint,
        },
        "data": data,
    }

    # Setup path for outputfile
    filepath = RAW_DIR / filename

    # Save meta data + API response
    with filepath.open("w", encoding="utf-8") as file:
        json.dump(
            payload,
            file,
            indent=2,
            ensure_ascii=False,
        )

    return filepath

def inspect(data: pd.DataFrame) -> pd.DataFrame:
    """Summarize the structure and completeness of the extracted data."""

    if data.empty:
        raise ValueError("K-Pro dataset is empty")

    print(f"No. of recods:    {len(data):,}")
    print(f"No. of columns: {len(data.columns):,}")
    
    summary = pd.DataFrame({
        "Column": data.columns,
        "Data type": data.dtypes.astype(str).values,
        "Non-null": data.notna().sum().values,
        "Missing %": np.round((data.isna().mean() * 100).values, 2),
        "Unique": data.nunique().values,
        "Duplicates": data.duplicated().sum(),
    })

    return summary
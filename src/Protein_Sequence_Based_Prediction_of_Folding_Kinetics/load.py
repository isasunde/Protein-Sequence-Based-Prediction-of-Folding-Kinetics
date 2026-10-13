import os
import pandas as pd
from mysql import connector
from mysql.connector import Error

from pathlib import Path
from dotenv import load_dotenv

PROJECT_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(PROJECT_ROOT / ".env")

QUERY_PATH = PROJECT_ROOT / "sql" / "analysis_dataset.sql"

TARGET = "ln_kf_h2o"

REQUIRED_COLUMNS = [
    "measurement_id",
    "experiment_id",
    "protein_id",
    "mutation_id",
    "reference_context",
    "mutation_context",
    TARGET,
]

def load_data(data, connection, table_name):
    """
    Load a cleaned DataFrame into an existing MySQL table.

    Parameters
    ----------
    data : pandas.DataFrame
        Cleaned dataset.
    connection : mysql.connector.MySQLConnection
        Active MySQL connection.
    table_name : str
        Destination table name.

    Returns
    -------
    int
        Number of rows inserted.
    """
    if data.empty:
        print("No records to load.")
        return 0

    # Convert pandas missing values (NaN, NaT, pd.NA) to Python None.
    data = data.astype(object).where(pd.notna(data), None)

    column_sql = ", ".join(f"`{col}`" for col in data.columns)
    placeholders = ", ".join(["%s"] * len(data.columns))

    query = (
        f"INSERT INTO `{table_name}` ({column_sql}) "
        f"VALUES ({placeholders})"
    )

    records = list(data.itertuples(index=False, name=None))
    cursor = connection.cursor()

    try:
        cursor.executemany(query, records)
        connection.commit()

        inserted = cursor.rowcount
        print(f"Successfully inserted {inserted} records into {table_name}.")
        return inserted

    except Error:
        connection.rollback()
        raise

    finally:
        cursor.close()

def fetch_dataset() -> pd.DataFrame:
    """Fetch the canonical analytical dataset from MySQL."""

    config = {
        "host": os.environ["MYSQL_HOST"],
        "port": int(os.getenv("MYSQL_PORT", "3306")),
        "user": os.environ["MYSQL_USER"],
        "password": os.environ["MYSQL_PASSWORD"],
        "database": os.getenv("MYSQL_DATABASE", "kpro"),
    }

    query = QUERY_PATH.read_text(encoding="utf-8")

    # The query file contains USE kpro;.
    # Remove it because the connection already selects the database.
    query = "\n".join(
        line for line in query.splitlines()
        if line.strip().upper() != "USE KPRO;"
    )

    connection = connector.connect(**config)

    try:
        return pd.read_sql_query(query, connection)
    finally:
        connection.close()


def validate_dataset(df: pd.DataFrame) -> dict:
    """Return a compact quality report for the analytical dataset."""

    missing_columns = sorted(set(REQUIRED_COLUMNS) - set(df.columns))

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    report = {
        "rows": len(df),
        "unique_measurements": df["measurement_id"].nunique(),
        "duplicate_measurement_ids": int(
            df["measurement_id"].duplicated().sum()
        ),
        "unique_proteins": df["protein_id"].nunique(),
        "unique_mutations": df["mutation_id"].nunique(),
        "missing_target": int(df[TARGET].isna().sum()),
        "missing_reference_context": int(
            df["reference_context"].isna().sum()
        ),
        "missing_mutation_context": int(
            df["mutation_context"].isna().sum()
        ),
    }

    if report["duplicate_measurement_ids"]:
        raise ValueError(
            "Dataset contains multiple rows per measurement_id."
        )

    return report
"""Clean, validate, and load K-Prot data into MySQL."""

import os
import logging
from pathlib import Path

import pandas as pd
import mysql.connector
from dotenv import load_dotenv

from Protein_Sequence_Based_Prediction_of_Folding_Kinetics import (
    prepare_columns,
    validate,
    clean,
    transform,
    load_data,
    add_sequences,
)

# Configuration
PROJECT_DIR = Path(__file__).resolve().parent
RAW_DIR = PROJECT_DIR / "data" / "raw"
SQL_DIR = PROJECT_DIR / "sql"
TABLE_NAME = "placeholder"

PROBLEM_COLUMNS = [
    "PFAM",
    "CATH",
    "BACKGROUND_MUT_PDB",
    "BACKGROUND_MUT_UNIPROT",
]

EXP_COLUMNS = [
    "PROTEIN", "SOURCE", "LENGTH", "UniProt", "PFAM", "CATH",
    "EC_NUMBER", "PDB_wild", "MUTATED_CHAIN", "MUTATION_UNIPROT",
    "MUTATION_PDB", "SEC_STR", "RSA", "T", "pH", "BUFFER_NAME",
    "BUFFER_CONC", "ION_NAME", "ADDITIVES", "MEASURE_KIN",
    "METHOD_KIN", "ln(kf)_H2O", "dln(kf)_H2O", "ln(ku)_H2O",
    "dln(ku)_H2O", "ln(kf)_DEN", "dln(kf)_DEN", "ln(ku)_DEN",
    "dln(ku)_DEN", "dGKIN_H2O", "ddGKIN_H2O", "mf", "mu",
    "CM_KIN", "beta-T", "PhiF_H2O", "STATE", "REVERSIBILITY",
    "AUTHOR_KIN", "REFERENCE_KIN", "YEAR_KIN", "PMID_KIN",
    "REVIEW_DATE", "BACKGROUND_MUT_UNIPROT", "BACKGROUND_MUT_PDB",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)


import json


def load_raw_data(raw_dir=RAW_DIR):
    """Load API responses and combine their data records."""

    files = sorted(raw_dir.glob("*.json"))

    if not files:
        raise FileNotFoundError(
            f"No JSON files found in {raw_dir}. "
            "Run 1_data_acquisition.py first."
        )

    frames = []

    for file in files:
        logger.info("Loading %s", file.name)

        try:
            with file.open("r", encoding="utf-8") as json_file:
                payload = json.load(json_file)

            # Check the structure created by save_raw_response().
            if not isinstance(payload, dict):
                raise ValueError("Expected a JSON object.")

            if "data" not in payload:
                raise ValueError(
                    "Missing 'data' field in saved API response."
                )

            response_data = payload["data"]

            if isinstance(response_data, list):
                records = response_data
            elif isinstance(response_data, dict):
                records = [response_data]
            else:
                raise ValueError(
                    "The 'data' field must contain a list or dictionary."
                )

            if not records:
                logger.warning("No records found in %s", file.name)
                continue

            if not all(isinstance(record, dict) for record in records):
                raise ValueError(
                    "Expected each data record to be a dictionary."
                )

            frame = pd.json_normalize(records)

            # Preserve source information for traceability.
            frame["_source_file"] = file.name

            metadata = payload.get("metadata", {})
            if metadata.get("endpoint"):
                frame["_source_endpoint"] = metadata["endpoint"]

            frames.append(frame)

        except (json.JSONDecodeError, OSError, ValueError) as exc:
            raise ValueError(
                f"Failed to load {file.name}: {exc}"
            ) from exc

    if not frames:
        raise ValueError(
            f"No usable records were found in {len(files)} JSON files."
        )

    data_full = pd.concat(frames, ignore_index=True)

    logger.info(
        "Loaded %d records from %d files.",
        len(data_full),
        len(frames),
    )

    return data_full


def prepare_raw_data(data_full):
    """Standardise columns and handle list-valued fields."""

    logger.info("Preparing columns...")

    data_full = prepare_columns(data_full)

    # Preserve missing values rather than converting them to "nan".
    for column in PROBLEM_COLUMNS:
        if column in data_full.columns:
            data_full[column] = data_full[column].map(
                lambda value: str(value)
                if isinstance(value, list)
                else value
            )

    # Retain the expected schema and column order.
    # Extra columns are excluded intentionally.
    data_full = data_full.reindex(columns=EXP_COLUMNS)

    return data_full


def run_etl(data_full):
    """Clean, enrich, transform, and validate the dataset."""

    logger.info("Cleaning data...")
    data_clean = clean(data_full)

    logger.info("Adding sequences...")
    data_clean = add_sequences(data_clean)

    logger.info("Transforming data...")
    data = transform(data_clean)

    logger.info("Validating data...")
    validation_result = validate(data)

    # If validate() returns False on failure, prevent loading.
    if validation_result is False:
        raise ValueError("Data validation failed; refusing database load.")

    logger.info("ETL complete: %d records ready.", len(data))

    return data

def get_connection():
    """Create a MySQL connection from environment variables."""

    load_dotenv(PROJECT_DIR / ".env")

    required = [
        "MYSQL_HOST",
        "MYSQL_PORT",
        "MYSQL_USER",
        "MYSQL_PASSWORD",
        "MYSQL_DATABASE",
    ]

    missing = [key for key in required if not os.getenv(key)]

    if missing:
        raise ValueError(
            f"Missing database configuration: {', '.join(missing)}"
        )

    return mysql.connector.connect(
        host=os.environ["MYSQL_HOST"],
        port=int(os.environ["MYSQL_PORT"]),
        user=os.environ["MYSQL_USER"],
        password=os.environ["MYSQL_PASSWORD"],
        database=os.environ["MYSQL_DATABASE"],
    )


def load_to_database(data):
    """Load validated data and apply database sorting."""

    connection = get_connection()

    try:
        logger.info("Loading %d records into MySQL...", len(data))

        load_data(
            data=data,
            connection=connection,
            table_name=TABLE_NAME,
        )

        # Ensure any pending changes are committed.
        connection.commit()

        sort_file = SQL_DIR / "sort_data.sql"

        if sort_file.exists():
            logger.info("Applying database sorting...")
            with sort_file.open(encoding="utf-8") as sql_file:
                with connection.cursor() as cursor:
                    cursor.execute(sql_file.read())

            connection.commit()
        else:
            logger.warning("Sorting SQL file not found: %s", sort_file)

        logger.info("Database load complete.")

    except Exception:
        connection.rollback()
        raise

    finally:
        connection.close()

def main():
    """Run the complete ETL pipeline."""

    data_full = load_raw_data()
    data_prepared = prepare_raw_data(data_full)
    data = run_etl(data_prepared)
    load_to_database(data)

if __name__ == "__main__":
    main()

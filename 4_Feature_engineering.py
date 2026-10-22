import os
import re
import numpy as np
import pandas as pd
import mysql.connector

from pathlib import Path
from getpass import getpass
from dotenv import load_dotenv
from Protein_Sequence_Based_Prediction_of_Folding_Kinetics import (
    fetch_dataset,
    data_embed,
    load_embeddings,
)

# ============================================
# Configuration
# ============================================
PROJECT_ROOT = Path(__file__).resolve().parent
load_dotenv(PROJECT_ROOT / ".env")

EMBEDDING_PATH = (
    PROJECT_ROOT / "data" / "processed" / "esm_embeddings.csv"
)

PERFORM_EMBEDDING = True
BATCH_SIZE = 8

TARGET = "ln_kf_h2o"
SEQUENCE = "sequence"

METADATA_COLUMNS = [
    "measurement_id",
    "experiment_id",
    "variable_id",
    "protein_id",
    "protein_name",
    "uniprot_id",
    "mutation_id",
    "mutation_uniprot",
    "mutation_pdb",
    "reference_context",
    "mutation_context",
]

# ============================================
# 1. Validate and prepare sequence data
# ============================================
def prepare_sequences(data: pd.DataFrame) -> pd.DataFrame:
    """Clean sequence contexts and validate the prediction target."""

    required = [
        "measurement_id",
        "protein_id",
        "mutation_id",
        "mutation_uniprot",
        'reference_context',
        'mutation_context',
        TARGET,
    ]

    missing = set(required) - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if data["measurement_id"].duplicated().any():
        raise ValueError(
            "SQL query returned duplicate measurement_id values. "
            "Check the normalized tables and joins before embedding."
        )

    data = data.copy()

    # Select the reference sequence for WT records and the mutant sequence otherwise.
    is_wildtype = (
        data["mutation_uniprot"]
        .astype("string")
        .str.fullmatch(r"\s*WT\s*", case=False, na=False)
    )

    data["sequence"] = data["mutation_context"].where(
        ~is_wildtype,
        data["reference_context"],
    )

    data["sequence"] = (
        data["sequence"]
        .astype("string")
        .str.upper()
        .str.replace(r"\s+", "", regex=True)
    )

    data[TARGET] = pd.to_numeric(data[TARGET], errors="coerce")

    valid_sequence = data['sequence'].str.fullmatch(
        r"[ACDEFGHIKLMNPQRSTVWYX]+",
        na=False,
    )

    valid_target = data[TARGET].notna()
    valid = valid_sequence & valid_target

    excluded = int((~valid).sum())
    if excluded:
        print(f"Excluded {excluded} rows with invalid sequences or missing targets.")

    data = data.loc[valid].reset_index(drop=True)

    if data.empty:
        raise ValueError("No valid observations remain after data preparation.")

    print(f"Valid observations: {len(data):,}")
    print(f"Unique proteins: {data['protein_id'].nunique():,}")
    print(f"Unique mutations: {data['mutation_id'].nunique():,}")
    print(f"Unique sequences: {data['sequence'].nunique():,}")

    return data


# ============================================
# 2. Generate or load embeddings
# ============================================
def load_or_create_embeddings(data: pd.DataFrame):
    """Generate embeddings or align a cached file to the current SQL dataset."""

    # Adapt SQL column names to the existing repository embedding interface.
    embedding_data = data.rename(
        columns={
            SEQUENCE: "sequence",
            TARGET: "ln(kf)_H2O",
        }
    )

    if PERFORM_EMBEDDING:
        token = getpass("Token from Forge: ")

        X, y, embeddings_df = data_embed(
            embedding_data,
            token,
            output_file=str(EMBEDDING_PATH),
            batch_size=BATCH_SIZE,
            sequence_col="sequence",
            target_col="ln(kf)_H2O",
            metadata_cols=METADATA_COLUMNS,
        )

    else:
        if not EMBEDDING_PATH.exists():
            raise FileNotFoundError(
                f"Cached embeddings not found: {EMBEDDING_PATH}\n"
                "Set PERFORM_EMBEDDING = True to generate them."
            )

        _, _, embeddings_df = load_embeddings(str(EMBEDDING_PATH))

        if "measurement_id" not in embeddings_df.columns:
            raise ValueError(
                "Cached embeddings lack measurement_id. "
                "Regenerate them using the updated data_embed() function."
            )

        embeddings_df["measurement_id"] = pd.to_numeric(
            embeddings_df["measurement_id"], errors="raise"
        ).astype("int64")

        if embeddings_df["measurement_id"].duplicated().any():
            raise ValueError("Cached embeddings contain duplicate measurement IDs.")

        expected_ids = data["measurement_id"].astype("int64").tolist()
        cached_ids = set(embeddings_df["measurement_id"])

        missing_ids = set(expected_ids) - cached_ids
        extra_ids = cached_ids - set(expected_ids)

        if missing_ids or extra_ids:
            raise ValueError(
                "Cached embeddings do not match the current SQL dataset. "
                f"Missing measurements: {len(missing_ids)}; "
                f"extra measurements: {len(extra_ids)}. "
                "Regenerate the embeddings."
            )

        # Reorder cached rows to match the current SQL data exactly.
        embeddings_df = (
            embeddings_df
            .set_index("measurement_id", drop=False)
            .loc[expected_ids]
            .reset_index(drop=True)
        )

        expected_sequences = data[SEQUENCE].tolist()

        if embeddings_df["sequence"].tolist() != expected_sequences:
            raise ValueError(
                "Cached sequences differ from SQL. Regenerate the embeddings."
            )

        if not np.allclose(
            embeddings_df["label"].to_numpy(dtype=float),
            data[TARGET].to_numpy(dtype=float),
            rtol=1e-5,
            atol=1e-6,
        ):
            raise ValueError(
                "Cached targets differ from SQL. Regenerate the embeddings."
            )

        feature_cols = [
            col for col in embeddings_df.columns
            if col.startswith("emb_") and col[4:].isdigit()
        ]

        X, y, _ = load_embeddings(str(EMBEDDING_PATH))

        # Rebuild X and y from the aligned dataframe to preserve row order.
        X = embeddings_df[feature_cols].to_numpy(dtype=np.float32)
        y = embeddings_df["label"].to_numpy(dtype=float)

    return X, y, embeddings_df


# ============================================
# 3. Run feature engineering
# ============================================
def main():
    data = fetch_dataset('feature_dataset.sql')
    data = prepare_sequences(data)

    X, y, embeddings_df = load_or_create_embeddings(data)

    if not (
        len(X) == len(y) == len(embeddings_df) == len(data)
    ):
        raise ValueError("Feature, target, and metadata lengths do not match.")

    print("\nFeature engineering complete.")
    print(f"Observations: {X.shape[0]:,}")
    print(f"Embedding dimensions: {X.shape[1]:,}")
    print(f"Target: {TARGET}")
    print(f"Saved embeddings: {EMBEDDING_PATH}")

    return X, y, embeddings_df


if __name__ == "__main__":
    X, y, embeddings_df = main()

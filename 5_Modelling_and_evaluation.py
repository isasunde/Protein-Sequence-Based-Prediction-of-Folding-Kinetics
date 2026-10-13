import os
from pathlib import Path

import mysql.connector
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dotenv import load_dotenv
from xgboost import XGBRegressor
from sklearn.linear_model import LinearRegression, Ridge

from Protein_Sequence_Based_Prediction_of_Folding_Kinetics import (
    load_embeddings,
    tsne_visualization,
    similarity_visualization,
    performance_histogram,
    true_vs_predicted_scatter,
    nested_cv,
    leave_one_protein_out_cv,
    evaluate_one_protein,
    summarize_cv_results,
)

# ============================================
# Configuration
# ============================================

PROJECT_ROOT = Path(__file__).resolve().parent
load_dotenv(PROJECT_ROOT / ".env")

EMBEDDING_PATH = (
    PROJECT_ROOT / "data" / "processed" / "esm_embeddings.csv"
)
RESULTS_DIR = PROJECT_ROOT / "reports" / "results"
FIGURES_DIR = PROJECT_ROOT / "reports" / "figures"

TARGET = "ln_kf_h2o"
RANDOM_STATE = 42

# Enable or disable computationally expensive steps.
PLOT_EMBEDDINGS = True
RUN_NESTED_CV = True
RUN_LOPO_CV = True
EVALUATE_ONE_PROTEIN = True

# ============================================
# 1. Retrieve metadata and target from SQL
# ============================================
def fetch_analysis_data() -> pd.DataFrame:
    """Retrieve measurement-level targets and metadata from MySQL."""

    config = {
        "host": os.environ["MYSQL_HOST"],
        "port": int(os.getenv("MYSQL_PORT", "3306")),
        "user": os.environ["MYSQL_USER"],
        "password": os.environ["MYSQL_PASSWORD"],
        "database": os.getenv("MYSQL_DATABASE", "kpro"),
    }

    query = """
        SELECT
            mv.variable_id,
            mv.measurement_id,
            mv.experiment_id,

            p.protein_id,
            p.protein_name,
            p.uniprot_id,

            mut.mutation_id,
            mut.mutation_uniprot,
            mut.mutation_pdb,

            mv.reference_context,
            mv.mutation_context,

            m.ln_kf_h2o,
            m.dln_kf_h2o,

            e.temperature,
            e.ph,
            e.buffer_name,
            e.buffer_conc,
            e.additives,
            e.kinetic_method,
            e.kinetic_measurement,
            e.state,
            e.reversibility,

            pub.publication_id,
            pub.pmid,
            pub.publication_year

        FROM model_variable AS mv

        JOIN measurement AS m
            ON m.measurement_id = mv.measurement_id
           AND m.experiment_id = mv.experiment_id

        JOIN experiment AS e
            ON e.experiment_id = mv.experiment_id
           AND e.protein_id = mv.protein_id
           AND e.mutation_id = mv.mutation_id

        JOIN protein AS p
            ON p.protein_id = mv.protein_id

        JOIN mutation AS mut
            ON mut.mutation_id = mv.mutation_id

        JOIN publication AS pub
            ON pub.publication_id = e.publication_id

        ORDER BY mv.measurement_id;
    """

    connection = mysql.connector.connect(**config)

    try:
        return pd.read_sql_query(query, connection)
    finally:
        connection.close()


# ============================================
# 2. Load embeddings and align with SQL data
# ============================================
def load_model_dataset():
    """
    Load cached ESM embeddings and join them with SQL metadata.

    The embedding file must contain measurement_id, sequence, label,
    and columns named emb_0, emb_1, ..., emb_n.
    """

    if not EMBEDDING_PATH.exists():
        raise FileNotFoundError(
            f"Embedding file not found: {EMBEDDING_PATH}\n"
            "Run 4_Feature_engineering.py first."
        )

    _, _, cached_embeddings = load_embeddings(str(EMBEDDING_PATH))
    sql_data = fetch_analysis_data()

    # Validate SQL data before merging.
    for name, frame in [
        ("SQL data", sql_data),
        ("Cached embeddings", cached_embeddings),
    ]:
        if "measurement_id" not in frame.columns:
            raise ValueError(f"{name} is missing measurement_id.")

        if frame["measurement_id"].isna().any():
            raise ValueError(f"{name} contains missing measurement IDs.")

        if frame["measurement_id"].duplicated().any():
            raise ValueError(
                f"{name} contains duplicate measurement IDs. "
                "Check the SQL joins or embedding generation."
            )

    # CSV files can infer integer IDs as floats; normalize both sides.
    sql_data["measurement_id"] = pd.to_numeric(
        sql_data["measurement_id"], errors="raise"
    ).astype("int64")

    cached_embeddings["measurement_id"] = pd.to_numeric(
        cached_embeddings["measurement_id"], errors="raise"
    ).astype("int64")

    # Identify feature columns explicitly. Metadata must not enter X.
    feature_cols = [
        col for col in cached_embeddings.columns
        if col.startswith("emb_") and col[4:].isdigit()
    ]

    if not feature_cols:
        raise ValueError("No ESM embedding columns were found.")

    required_cache_cols = {"sequence", "label", "measurement_id"}
    missing = required_cache_cols - set(cached_embeddings.columns)

    if missing:
        raise ValueError(
            f"Cached embeddings are missing columns: {sorted(missing)}. "
            "Regenerate them with 4_Feature_engineering.py."
        )

    # Keep only the embeddings and fields needed to validate alignment.
    embedding_data = cached_embeddings[
        ["measurement_id", "sequence", "label", *feature_cols]
    ].copy()

    # SQL is the source of truth for targets and metadata.
    data = sql_data.merge(
        embedding_data,
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    )

    if data[TARGET].isna().any():
        raise ValueError("SQL dataset contains missing target values.")

    if data["sequence"].isna().any() or data["mutation_context"].isna().any():
        raise ValueError("Sequence data contain missing values.")

    # Check that the cached sequence matches the SQL sequence context.
    sql_sequences = (
        data["mutation_context"].astype(str).str.upper().str.replace(
            r"\s+", "", regex=True
        )
    )
    cached_sequences = (
        data["sequence"].astype(str).str.upper().str.replace(
            r"\s+", "", regex=True
        )
    )

    if not sql_sequences.equals(cached_sequences):
        raise ValueError(
            "Cached sequences do not match mutation_context in SQL. "
            "Regenerate the embeddings."
        )

    # Ensure the cached labels correspond to the current SQL targets.
    labels = pd.to_numeric(data["label"], errors="coerce").to_numpy()
    targets = pd.to_numeric(data[TARGET], errors="coerce").to_numpy()

    if not np.isfinite(labels).all() or not np.allclose(
        labels, targets, rtol=1e-5, atol=1e-6
    ):
        raise ValueError(
            "Cached labels do not match SQL targets. "
            "Regenerate the embeddings."
        )

    # Use protein name for grouping
    data["protein"] = data["protein_name"].astype(str)

    # Retain a useful mutation label for the existing visualizations.
    data["mutation"] = (
        data["mutation_pdb"]
        .fillna(data["mutation_uniprot"])
        .fillna("Unknown")
    )

    X = data[feature_cols].to_numpy(dtype=np.float32)
    y = targets.astype(float)

    if not np.isfinite(X).all():
        raise ValueError("Embedding matrix contains NaN or infinite values.")

    print("Dataset loaded successfully.")
    print(f"Measurements: {len(data):,}")
    print(f"Unique proteins: {data['protein'].nunique():,}")
    print(f"Unique mutations: {data['mutation_id'].nunique():,}")
    print(f"Embedding dimensions: {X.shape[1]:,}")

    return X, y, data


# ============================================
# 3. Define regression models
# ============================================

def define_models():
    """Return the regression models and their scaling configurations."""

    return {
        "LinearRegression": {
            "model": LinearRegression(),
            "scale": True,
        },
        "Ridge": {
            "model": Ridge(alpha=10),
            "scale": True,
        },
        "XGBoost": {
            "model": XGBRegressor(
                n_estimators=250,
                max_depth=4,
                learning_rate=0.05,
                subsample=0.9,
                colsample_bytree=0.9,
                random_state=RANDOM_STATE,
                n_jobs=-1,
            ),
            "scale": False,
        },
    }


# ============================================
# 4. Explore embeddings
# ============================================

def explore_embeddings(X, y, data):
    """Visualize embedding structure and sequence similarity."""

    if not PLOT_EMBEDDINGS:
        return

    fig, ax = tsne_visualization(X, y)
    fig.savefig(
        FIGURES_DIR / "esm_tsne.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()

    fig, ax = similarity_visualization(X, y, data)
    fig.savefig(
        FIGURES_DIR / "esm_similarity.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


# ============================================
# 5. Hyperparameter tuning
# ============================================

def tune_models(X, y, data, models):
    """Run nested cross-validation for Ridge and XGBoost."""

    tuning_results = {}

    tuning_config = {
        "XGBoost": {
            "parameter": "n_estimators",
            "values": [100, 250, 500, 750, 1000],
        },
        "Ridge": {
            "parameter": "alpha",
            "values": [0.01, 0.1, 1, 10, 100],
        },
    }

    for model_name, config in tuning_config.items():
        print(f"\nTuning {model_name}...")

        outer_results, inner_results = nested_cv(
            X,
            y,
            data,
            model_config=models[model_name],
            hyperparameter_name=config["parameter"],
            hyperparameters=config["values"],
        )

        outer_results.to_csv(
            RESULTS_DIR / f"{model_name.lower()}_nested_cv.csv",
            index=False,
        )
        inner_results.to_csv(
            RESULTS_DIR / f"{model_name.lower()}_inner_cv.csv",
            index=False,
        )

        tuning_results[model_name] = {
            "outer": outer_results,
            "inner": inner_results,
        }

    return tuning_results


# ============================================
# 6. Compare models using leave-one-protein-out CV
# ============================================

def compare_models(X, y, data, models):
    """Evaluate model generalization to held-out proteins."""

    print("\nRunning leave-one-protein-out cross-validation...")

    results = leave_one_protein_out_cv(
        X,
        y,
        data,
        models,
    )

    summary = summarize_cv_results(results, models)

    results.to_csv(
        RESULTS_DIR / "leave_one_protein_out_results.csv",
        index=False,
    )
    summary.to_csv(
        RESULTS_DIR / "model_performance_summary.csv",
        index=False,
    )

    print("\nModel performance summary:")
    print(summary.to_string(index=False))

    fig, ax = performance_histogram(results)
    fig.savefig(
        FIGURES_DIR / "model_performance.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()

    return results, summary


# ============================================
# 7. Evaluate an individual protein
# ============================================

def evaluate_protein(X, y, data, models):
    """Visualize predictions for a selected held-out protein."""

    # Select a protein with the most observations, making the evaluation
    # plot more informative than choosing an arbitrary list index.
    protein_counts = data["protein"].value_counts()

    if protein_counts.empty:
        print("No proteins available for individual evaluation.")
        return

    protein = protein_counts.index[0]
    protein_name = data.loc[
        data["protein"] == protein, "protein_name"
    ].dropna()

    display_name = (
        f"{protein_name.iloc[0]} (ID: {protein})"
        if not protein_name.empty
        else protein
    )

    print(
        f"\nEvaluating held-out protein: {display_name} "
        f"({protein_counts.loc[protein]} measurements)"
    )

    y_test, y_pred, mse = evaluate_one_protein(
        data,
        X,
        y,
        models,
        protein,
    )

    fig, axes = true_vs_predicted_scatter(
        y_test,
        y_pred,
        mse,
        models,
        display_name,
    )

    fig.savefig(
        FIGURES_DIR / "individual_protein_predictions.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


# ============================================
# 8. Run pipeline
# ============================================

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    X, y, data = load_model_dataset()
    models = define_models()

    explore_embeddings(X, y, data)

    if RUN_NESTED_CV:
        tune_models(X, y, data, models)

    if RUN_LOPO_CV:
        compare_models(X, y, data, models)

    if EVALUATE_ONE_PROTEIN:
        evaluate_protein(X, y, data, models)

    print("\nModelling and evaluation complete.")


if __name__ == "__main__":
    main()

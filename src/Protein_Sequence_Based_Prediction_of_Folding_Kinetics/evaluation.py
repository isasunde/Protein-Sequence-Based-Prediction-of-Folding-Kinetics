import numpy as np
import pandas as pd

from sklearn.base import clone
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

def leave_one_protein_out_cv(
    X,
    y,
    embeddings_df,
    models,
):
    """
    Perform leave-one-protein-out cross-validation.

    Each unique protein is left out once as the test set, while all
    remaining proteins are used for training.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix.

    y : np.ndarray
        Target values.

    embeddings_df : pd.DataFrame
        DataFrame containing a 'protein' column corresponding to the
        samples in X and y.

    models : dict
        Dictionary containing model configurations. Each model must
        contain a 'model' and 'scale' entry.

    Returns
    -------
    results_df : pd.DataFrame
        Fold-level performance metrics for each model.
    """

    proteins = embeddings_df["protein"].unique()
    results = []

    for fold, protein in enumerate(
        tqdm(proteins, desc="Leave-one-protein-out CV"),
        start=1,
    ):
        test_mask = embeddings_df["protein"] == protein

        X_train = X[~test_mask]
        X_test = X[test_mask]

        y_train = y[~test_mask]
        y_test = y[test_mask]

        fold_results = {
            "Fold": fold,
            "Protein": protein,
        }

        for model_name, model_config in models.items():

            model = clone(model_config["model"])

            X_train_model = X_train
            X_test_model = X_test

            if model_config["scale"]:
                scaler = StandardScaler()

                X_train_model = scaler.fit_transform(
                    X_train
                )

                X_test_model = scaler.transform(
                    X_test
                )

            model.fit(
                X_train_model,
                y_train,
            )

            y_pred = model.predict(
                X_test_model
            )

            mse = mean_squared_error(
                y_test,
                y_pred,
            )

            rmse = np.sqrt(mse)

            r2 = r2_score(
                y_test,
                y_pred,
            )

            fold_results[f"{model_name}_MSE"] = mse
            fold_results[f"{model_name}_RMSE"] = rmse
            fold_results[f"{model_name}_R2"] = r2

        results.append(fold_results)

    results_df = pd.DataFrame(results)

    return results_df


def summarize_cv_results(results_df, models):
    """
    Calculate mean cross-validation performance for each model.

    Parameters
    ----------
    results_df : pd.DataFrame
        Fold-level cross-validation results.

    models : dict
        Dictionary of model configurations.

    Returns
    -------
    summary_df : pd.DataFrame
        Mean performance across all proteins.
    """

    summary = []

    for model_name in models:
        summary.append(
            {
                "Model": model_name,
                "MSE": results_df[
                    f"{model_name}_MSE"
                ].mean(),
                "RMSE": results_df[
                    f"{model_name}_RMSE"
                ].mean(),
                "R2": results_df[
                    f"{model_name}_R2"
                ].mean(),
            }
        )

    return pd.DataFrame(summary)

def nested_cv(
    X,
    y,
    embeddings_df,
    model_config,
    hyperparameter_name,
    hyperparameters,
    n_outer_splits=5,
    n_inner_splits=3,
    random_state=42,
):
    """
    Perform nested cross-validation with protein-level splitting.

    The outer loop estimates generalisation performance on unseen
    proteins, while the inner loop selects the optimal value of a
    model hyperparameter.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix with shape (n_samples, n_features).

    y : np.ndarray
        Target values.

    embeddings_df : pd.DataFrame
        DataFrame containing supplementary information about the
        embeddings. Must contain a 'protein' column corresponding
        to the rows in X and y.

    model_config : dict
        Model configuration containing a model instance and a
        scaling setting. For example:

        {
            "model": XGBRegressor(...),
            "scale": False,
        }

    hyperparameter_name : str
        Name of the hyperparameter being optimised.

    hyperparameters : list
        Values to test for the selected hyperparameter.

    n_outer_splits : int, optional
        Number of outer cross-validation folds. Default is 5.

    n_inner_splits : int, optional
        Number of inner cross-validation folds. Default is 3.

    random_state : int, optional
        Random seed for reproducibility. Default is 42.

    Returns
    -------
    results_df : pd.DataFrame
        Outer-fold performance and selected hyperparameter values.

    inner_results_df : pd.DataFrame
        Inner cross-validation results for all hyperparameter values.
    """

    proteins = embeddings_df["protein"].to_numpy()

    unique_proteins = np.unique(proteins)

    if n_outer_splits > len(unique_proteins):
        raise ValueError(
            "n_outer_splits cannot be greater than the number "
            "of unique proteins."
        )

    results = []
    inner_results = []

    # ------------------------------------------------------------
    # Create reproducible outer protein-level folds
    # ------------------------------------------------------------
    rng = np.random.default_rng(random_state)
    shuffled_proteins = rng.permutation(unique_proteins)

    outer_fold_sizes = np.full(
        n_outer_splits,
        len(shuffled_proteins) // n_outer_splits,
    )

    outer_fold_sizes[
        :len(shuffled_proteins) % n_outer_splits
    ] += 1

    outer_folds = []
    start = 0

    for size in outer_fold_sizes:
        outer_folds.append(
            shuffled_proteins[start:start + size]
        )
        start += size

    # ------------------------------------------------------------
    # Outer CV
    # ------------------------------------------------------------
    for fold, test_proteins in enumerate(
        outer_folds,
        start=1,
    ):
        print("\n==============================")
        print(
            f"Outer Fold {fold}/{n_outer_splits}"
        )
        print(
            f"Test proteins: {list(test_proteins)}"
        )
        print("==============================\n")

        # Outer test split
        test_mask = np.isin(
            proteins,
            test_proteins,
        )

        X_test = X[test_mask]
        y_test = y[test_mask]

        X_train_val = X[~test_mask]
        y_train_val = y[~test_mask]

        groups_train_val = proteins[~test_mask]

        # --------------------------------------------------------
        # Inner CV
        # --------------------------------------------------------
        inner_cv = GroupKFold(
            n_splits=n_inner_splits,
        )

        parameter_mse = {
            parameter: []
            for parameter in hyperparameters
        }

        for parameter in tqdm(
            hyperparameters,
            desc="Hyperparameter tuning",
        ):
            for inner_train_idx, inner_val_idx in inner_cv.split(
                X_train_val,
                y_train_val,
                groups=groups_train_val,
            ):
                X_train = X_train_val[inner_train_idx]
                X_val = X_train_val[inner_val_idx]

                y_train = y_train_val[inner_train_idx]
                y_val = y_train_val[inner_val_idx]

                # ------------------------------------------------
                # Optional feature scaling
                # ------------------------------------------------
                if model_config["scale"]:
                    scaler = StandardScaler()

                    X_train = scaler.fit_transform(
                        X_train
                    )

                    X_val = scaler.transform(
                        X_val
                    )

                # ------------------------------------------------
                # Configure model
                # ------------------------------------------------
                model = clone(
                    model_config["model"]
                )

                model.set_params(
                    **{
                        hyperparameter_name: parameter
                    }
                )

                # ------------------------------------------------
                # Train and evaluate
                # ------------------------------------------------
                model.fit(
                    X_train,
                    y_train,
                )

                predictions = model.predict(
                    X_val
                )

                mse = mean_squared_error(
                    y_val,
                    predictions,
                )

                parameter_mse[parameter].append(
                    mse
                )

                inner_results.append(
                    {
                        "Outer_Fold": fold,
                        hyperparameter_name: parameter,
                        "MSE": mse,
                    }
                )

        # --------------------------------------------------------
        # Select best hyperparameter
        # --------------------------------------------------------
        average_mse = {
            parameter: np.mean(mse_values)
            for parameter, mse_values
            in parameter_mse.items()
        }

        best_parameter = min(
            average_mse,
            key=average_mse.get,
        )

        print(
            "Inner CV validation MSE:"
        )

        for parameter in hyperparameters:
            print(
                f"  {hyperparameter_name}={parameter}: "
                f"{average_mse[parameter]:.4f}"
            )

        print(
            f"\nBest {hyperparameter_name}: "
            f"{best_parameter}"
        )

        # --------------------------------------------------------
        # Retrain using complete outer training set
        # --------------------------------------------------------
        X_train_final = X_train_val.copy()
        X_test_final = X_test.copy()

        if model_config["scale"]:
            scaler = StandardScaler()

            X_train_final = scaler.fit_transform(
                X_train_final
            )

            X_test_final = scaler.transform(
                X_test_final
            )

        final_model = clone(
            model_config["model"]
        )

        final_model.set_params(
            **{
                hyperparameter_name: best_parameter
            }
        )

        final_model.fit(
            X_train_final,
            y_train_val,
        )

        predictions_test = final_model.predict(
            X_test_final
        )

        test_mse = mean_squared_error(
            y_test,
            predictions_test,
        )

        results.append(
            {
                "Fold": fold,
                "Test_Proteins": list(test_proteins),
                "MSE": test_mse,
                hyperparameter_name: best_parameter,
            }
        )

        print(
            f"Outer fold {fold} test MSE: "
            f"{test_mse:.4f}"
        )

    # ------------------------------------------------------------
    # Convert results to DataFrames
    # ------------------------------------------------------------
    results_df = pd.DataFrame(results)
    inner_results_df = pd.DataFrame(inner_results)

    print("\n==============================")
    print("Final Nested CV Results")
    print("==============================\n")

    print(
        results_df.to_string(index=False)
    )

    mean_test_mse = results_df["MSE"].mean()

    print(
        f"\nOverall Nested CV Generalisation MSE: "
        f"{mean_test_mse:.4f}"
    )

    return results_df, inner_results_df

def evaluate_one_protein(
    embeddings_df: pd.DataFrame,
    X: np.ndarray,
    y: np.ndarray,
    models: dict,
    protein: str,
):
    """
    Train and evaluate models using one protein as the test set.

    Parameters
    ----------
    embeddings_df : pd.DataFrame
        DataFrame containing protein information for each sample.
    X : np.ndarray
        Feature matrix.
    y : np.ndarray
        Target values.
    models : dict
        Dictionary containing model configurations. Each configuration
        should contain a ``model`` and ``scale`` entry.
    protein : str
        Protein to use as the test set.

    Returns
    -------
    y_test : np.ndarray
        True target values for the test protein.
    y_pred : dict
        Predicted target values for each model.
    mse : dict
        Mean squared error for each model.
    """

    # Split data into training and test sets
    test_mask = embeddings_df["protein"] == protein

    X_train = X[~test_mask]
    X_test = X[test_mask]
    y_train = y[~test_mask]
    y_test = y[test_mask]

    # Set up result storage
    y_pred = {}
    mse = {}

    # Train and evaluate each model
    for model_name, model_config in models.items():
        model = clone(model_config["model"])

        if model_config["scale"]:
            scaler = StandardScaler()
            X_train_model = scaler.fit_transform(X_train)
            X_test_model = scaler.transform(X_test)
        else:
            X_train_model = X_train
            X_test_model = X_test

        model.fit(X_train_model, y_train)

        y_pred[model_name] = model.predict(X_test_model)

        mse[model_name] = mean_squared_error(
            y_test,
            y_pred[model_name],
        )

        print(
            f"{model_name} MSE on protein {protein}: "
            f"{mse[model_name]:.4f}"
        )

    return y_test, y_pred, mse
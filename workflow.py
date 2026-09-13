# Housekeeping
import pandas as pd
import matplotlib.pyplot as plt

from getpass import getpass
from xgboost import XGBRegressor
from sklearn.linear_model import LinearRegression, Ridge
from folding_kinetics import (
    data_filter,
    data_embed,
    load_embeddings,
    tsne_visualization,
    similarity_visualization,
    performance_histogram,
    true_vs_predicted_scatter,
    nested_cv,
    leave_one_protein_out_cv,
    evaluate_one_protein,
    summarize_cv_results
) 

# Inputs
file_path = "data/kpro.tsv" # path to the dataset
perform_embedding = False # whether to perform embedding or load the saved embeddings

# ============================================
# Load dataset and preprocess
# ============================================
# Load the data set 
data = pd.read_csv(file_path, sep='\t')

# Attributes to be included
columns_of_interest = ['sequence','sequence_mutated', 'ln(kf)_H2O']

# Filter data 
data = data_filter(data, columns_of_interest)

# Display properties of dataset
print('The size of the dataset is:', data.shape)
print('The number of attributes is:', data.shape[1])
print('The number of samples is:', data.shape[0])

# ============================================
# Feature Engineering - ESM representations
# ============================================
if perform_embedding:
    # Embed sequences using ESM3
    token = getpass("Token from Forge: ") # token to Hugging Face to access weights of the pretrained ESM models
    X, y, embeddings_df = data_embed(data, token, output_file="data/esm_embeddings.csv", batch_size=8)

else:
    X, y, embeddings_df = load_embeddings(
        "data/esm_embeddings_final.csv"
    )

    data = pd.read_csv(file_path, sep="\t")

    kpro_info = (
        data[
            ['sequence_mutated', 'PROTEIN', 'MUTATION_PDB']
        ]
        .drop_duplicates(subset='sequence_mutated')
        .rename(
            columns={
                'sequence_mutated': 'sequence',
                'PROTEIN': 'protein',
                'MUTATION_PDB': 'mutation',
            }
        )
    )

    embeddings_df = embeddings_df.merge(
        kpro_info[['sequence', 'protein', 'mutation']],
        on='sequence',
        how='left',
        validate='many_to_one',
    )

# Dimensionality of data
N, M = X.shape

# Define dimension for ANN
input_dim = M
output_dim = 1

print(f"X samples: {len(X)}")
print(f"y samples: {len(y)}")
print(f"Embedding metadata rows: {len(embeddings_df)}")

# ============================================
# Explore embeddings
# ============================================
# Display t-SNE of ESM embeddings
fig, ax = tsne_visualization(X, y)
plt.show()

# Display similarity matrix of ESM embeddings
fig, ax = similarity_visualization(X, y, embeddings_df)
plt.show()

# ============================================
# Define models
# ============================================
models = {
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
            random_state=42,
            n_jobs=-1,
        ),
        "scale": False,
    },
}

# ============================================
# Model selection and hyperparameter tuning
# ============================================
# Tune XGBoost
xgb_results, xgb_inner_results = nested_cv(
    X,
    y,
    embeddings_df,
    model_config=models["XGBoost"],
    hyperparameter_name="n_estimators",
    hyperparameters=[100, 250, 500, 750, 1000],
)

# Tune Ridge
ridge_results, ridge_inner_results = nested_cv(
    X,
    y,
    embeddings_df,
    model_config=models["Ridge"],
    hyperparameter_name="alpha",
    hyperparameters=[0.01, 0.1, 1, 10, 100],
)

# ============================================
# Compare model performance
# ============================================
results_df = leave_one_protein_out_cv(
    X,
    y,
    embeddings_df,
    models,
)

summary_df = summarize_cv_results(results_df, models)

# Plot model performance
fig, ax = performance_histogram(results_df)
plt.show()

# ============================================
# Evaluate model on a specific protein
# ============================================
# Choose a protein to visualize
protein = embeddings_df["protein"].unique()[8]

y_test, y_pred, mse = evaluate_one_protein(
    embeddings_df,
    X,
    y,
    models,
    protein,
)

fig, axs = true_vs_predicted_scatter(
    y_test,
    y_pred,
    mse,
    models,
    protein,
)
plt.show()

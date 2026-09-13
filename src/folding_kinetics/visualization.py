import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.manifold import TSNE

def tsne_visualization(
        X: np.ndarray,
        y: np.ndarray, 
        perplexity: float=40.0, 
        random_state: int=1234
    ):
    """
    Visualizes the embeddings using t-SNE.

    Parameters:
    - X: The embeddings (numpy array).
    - y: The target values (numpy array).
    """
    
    # If using a small dataset, set perplexity to N/5 as it should be less than N/3 (should be in range 5-50)
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=random_state)
    X_tsne = tsne.fit_transform(X)

    fig, ax = plt.subplots(
            figsize=(10, 6)
        )
    sc = ax.scatter(
        X_tsne[:, 0],
        X_tsne[:, 1],
        c=y,
        cmap="viridis",
        alpha=0.7,
        s=12
    )

    fig.colorbar(sc, ax=ax, label="y value")
    ax.set_title("t-SNE of ESM embeddings")
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")

    return fig, ax

def similarity_visualization(
        X: np.ndarray, 
        y: np.ndarray, 
        df: pd.DataFrame
    ):
    """
    Visualizes the similarity matrix using a heatmap.

    Parameters:
    - X: The embeddings (numpy array).
    - y: The target values (numpy array).
    - df: Supplementary information (pandas dataframe).
    """
    
    # Compute similarity matrix using cosine similarity
    similarity_matrix = np.dot(X, X.T)/(np.linalg.norm(X, axis=1)[:, None] * np.linalg.norm(X, axis=1)[None, :])
    #similarity_matrix = np.log(similarity_matrix)

    # Plot the similarity matrix
    fig, ax = plt.subplots(
                figsize=(10, 6)
            )
    
    sm = ax.imshow(similarity_matrix, cmap="Blues", interpolation="None")
    fig.colorbar(sm, ax=ax)
    #plt.xticks(embeddings_df['protein'].unique())
    #plt.yticks(embeddings_df['protein'].unique())
    ax.set_xlabel("Proteins")
    ax.set_ylabel("Proteins")
    ax.set_title("Similarity matrix of ESM representations")
    ax.grid(False)

    # Determine proteins with the lowest general similarity to others
    avg_similarity = np.mean(similarity_matrix, axis=1)
    protein_similarity = pd.DataFrame({
        'Protein': df['protein'],
        'Avg_Similarity': avg_similarity
    })
    protein_similarity = protein_similarity.drop_duplicates().sort_values(by='Avg_Similarity')
    print("\nProteins with lowest average similarity to others:")
    print(protein_similarity.head(10))
    print("\nFirst 10 unique protein names:")
    for protein in protein_similarity['Protein'].unique()[0:10]:
        print(protein)

    return fig, ax

def performance_histogram(
                          results_df: pd.DataFrame, 
                          ):
    """
    Visualizes the distribution of Mean Squared Error (MSE) across proteins for Ridge Regression and XGBoost models."""

    # Vizualize the results in a histogram
    fig, axs = plt.subplots(1,2, figsize=(12, 5))
    fig.suptitle(f'Distribution of MSE across Proteins')
    axs[0].hist(results_df['Ridge_MSE'], alpha=0.5, label='Ridge Regression MSE', bins=10)
    axs[0].set_xlabel("Mean Squared Error (MSE)")
    axs[0].set_ylabel("Frequency")
    axs[0].set_title(f'Ridge Regression (MSE: {np.mean(results_df['Ridge_MSE']):.2f})')

    axs[1].hist(results_df['XGBoost_MSE'], alpha=0.5, label='XGBoost MSE', bins=10)
    axs[1].set_xlabel("Mean Squared Error (MSE)")
    axs[1].set_ylabel("Frequency")
    axs[1].set_title(f'XGBoost (MSE: {np.mean(results_df['XGBoost_MSE']):.2f})')
    return fig, axs

def true_vs_predicted_scatter(
        y_test: np.ndarray,
        y_pred: dict,
        mse: dict, 
        models: dict,
        protein: str
    ):
    """
    Visualizes the true vs predicted values for a given model.

    Parameters:
    - results_df: DataFrame containing true and predicted values.
    - model_name: Name of the model to visualize ('Ridge' or 'XGBoost').
    """

    # Plot true vs predicted values
    fig, axs = plt.subplots(1,len(models), figsize=(len(models)*6, 5))
    fig.suptitle(f'Predicted vs. true on protein: {protein}')

    for i, model_name, ___ in enumerate(models.items()):
        axs[i].scatter(y_test, y_pred[model_name], alpha=0.7)
        axs[i].plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--')
        axs[i].set_xlabel("True values")
        axs[i].set_ylabel("Predicted values")
        axs[i].set_title(f'{model_name} (MSE: {mse[model_name]:.2f})')
    
    return fig, axs
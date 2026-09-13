'Functions required to embed sequences with ESM3 (COPIED from esm github tutorial)'
import torch
import pandas as pd

from tqdm import tqdm
from esm.sdk import client
from typing import Sequence
from concurrent.futures import ThreadPoolExecutor
from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    ESMProteinError,
    LogitsConfig,
    LogitsOutput,
    ProteinType,
)

EMBEDDING_CONFIG = LogitsConfig(
    sequence=True, return_embeddings=True, return_hidden_states=True
)

def embed_sequence(model: ESM3InferenceClient, sequence: str) -> LogitsOutput:
    protein = ESMProtein(sequence=sequence)
    protein_tensor = model.encode(protein)
    output = model.logits(protein_tensor, EMBEDDING_CONFIG)
    return output

def batch_embed(
    model: ESM3InferenceClient, inputs: Sequence[ProteinType]
) -> Sequence[LogitsOutput]:
    """Forge supports auto-batching. So batch_embed() is as simple as running a collection    
    of embed calls in parallel using asyncio.
    """
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(embed_sequence, model, protein) for protein in inputs
        ]
        results = []
        for future in futures:
            try:
                results.append(future.result())
            except Exception as e:
                results.append(ESMProteinError(500, str(e)))
    return results

def data_embed(
        data: pd.DataFrame, 
        token: str, 
        output_file: str="", 
        batch_size: int=8
    ):
    
    print("Embedding sequences using ESM3...")
    # Load pretrained ESM model
    model = client(
        model="esmc-300m-2024-12", url="https://forge.evolutionaryscale.ai", token=token
    )
    
    # Embed sequences - Using the mean pooled embeddings from the second to last layer
    NUM_ESMC_6B_LAYERS = 30-1   # ESMC-300m-2024-12 has 30 layers
    ONE_LAYER_EMBEDDING_CONFIG = LogitsConfig(
        return_hidden_states=True,
        ith_hidden_layer=NUM_ESMC_6B_LAYERS
    )
    #MEAN_POOLED_EMBEDDING_CONFIG = LogitsConfig(return_hidden_states=True)
    
    sequences = data["sequence_mutated"][0:10].tolist()
    
    all_mean_embeddings = []
    
    for i in tqdm(range(0, len(sequences), batch_size), desc="Embedding proteins"):
        batch = sequences[i:i + batch_size]
    
        outputs = batch_embed(
            model,
            batch
        )
    
        batch_embeddings = [
            torch.mean(output.hidden_states.float(), dim=-2).squeeze()
            for output in outputs
        ]
    
        all_mean_embeddings.extend(batch_embeddings)
    
    # ============================================
    # Prepare the data embeddings and targets
    # ============================================
    X = (
        torch.stack([embedding[-1, :] for embedding in all_mean_embeddings])
        .float()
        .numpy()
        )
    y = data['ln(kf)_H2O'][0:10].values

    if output_file:
        # Save embeddings and targets as CSV
        df = pd.DataFrame(X, columns=[f"emb_{i}" for i in range(X.shape[1])])
        df["label"] = y
        df["sequence"] = data["sequence_mutated"][0:len(df)].values
        df.to_csv(output_file, index=False)

    return X, y, df

def load_embeddings(file_path: str):
    """
    Load embeddings and targets from a CSV file.

    Parameters:
    file_path (str): Path to the CSV file containing embeddings and targets.

    Returns:
    tuple: A tuple containing the embeddings (X) and targets (y).
    """
    
    df = pd.read_csv(file_path)
    X = df.drop(columns=["label", "sequence"]).values
    y = df["label"].values
    return X, y, df
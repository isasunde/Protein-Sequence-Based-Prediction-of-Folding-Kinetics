'Functions required to embed sequences with ESM3 embed_sequence, batch_embed and _pool_embeding from esm github tutorial'

import numpy as np
import pandas as pd
import torch

from pathlib import Path
from tqdm import tqdm
from esm.sdk import client
from esm.sdk.api import (
    ESM3InferenceClient,
    ESMProtein,
    ESMProteinError,
    LogitsConfig,
    LogitsOutput,
)
from esm.sdk.forge import ESMCForgeInferenceClient
from concurrent.futures import ThreadPoolExecutor
from typing import Sequence

EMBEDDING_CONFIG = LogitsConfig(
    sequence=True, return_embeddings=True, return_hidden_states=True
)

def embed_sequence(
    model: ESM3InferenceClient,
    sequence: str,
    config: LogitsConfig,
) -> LogitsOutput:
    """Generate an ESM sequence representation."""
    protein = ESMProtein(sequence=sequence)
    protein_tensor = model.encode(protein)
    return model.logits(protein_tensor, config)

def batch_embed(
    model: ESM3InferenceClient,
    inputs: Sequence[str],
    config: LogitsConfig,
) -> Sequence[LogitsOutput]:
    """Embed a batch of sequences."""
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(embed_sequence, model, seq, config)
            for seq in inputs
        ]
        results = []

        for future in futures:
            try:
                results.append(future.result())
            except Exception as exc:
                results.append(ESMProteinError(500, str(exc)))

    return results

def _pool_embedding(output: LogitsOutput) -> np.ndarray:
    """Mean-pool hidden states across the sequence dimension."""
    if isinstance(output, ESMProteinError):
        raise RuntimeError(f"ESM embedding failed: {output}")

    hidden = output.hidden_states

    if hidden is None:
        raise ValueError("ESM returned no hidden states.")

    hidden = hidden.float()

    # Handle a possible layer dimension, followed by batch/sequence dimensions.
    if hidden.ndim == 4:
        hidden = hidden[-1]

    if hidden.ndim == 3:
        if hidden.shape[0] == 1:
            hidden = hidden.squeeze(0)
        else:
            raise ValueError(
                f"Unexpected hidden-state shape: {tuple(hidden.shape)}"
            )

    if hidden.ndim != 2:
        raise ValueError(
            f"Expected sequence x embedding dimensions; got {tuple(hidden.shape)}"
        )

    # Mean across residues; output is one vector per sequence.
    return hidden.mean(dim=0).cpu().numpy()


def data_embed(
    data: pd.DataFrame,
    token: str,
    output_file: str = "",
    batch_size: int = 8,
    sequence_col: str = "sequence_mutated",
    target_col: str = "ln(kf)_H2O",
    metadata_cols: list[str] | None = None,
):
    """
    Embed unique sequences and map their embeddings back to all observations.

    The output contains embedding columns, target, sequence, and selected metadata.
    """
    required = [sequence_col, target_col]
    missing = set(required) - set(data.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data = data.reset_index(drop=True).copy()

    if data[sequence_col].isna().any() or data[target_col].isna().any():
        raise ValueError("Remove rows with missing sequences or targets first.")

    # Preserve the row-to-measurement mapping, but avoid redundant ESM calls.
    sequences = data[sequence_col].astype(str).tolist()
    unique_sequences = list(dict.fromkeys(sequences))

    model = ESMCForgeInferenceClient(
        model="esmc-6b-2024-12",
        url="https://biohub.ai",
        token=token,
    )

    # Select the penultimate hidden layer.
    config = LogitsConfig(
        return_hidden_states=True,
        ith_hidden_layer=29,
    )

    sequence_embeddings = {}

    for start in tqdm(
        range(0, len(unique_sequences), batch_size),
        desc="Embedding unique sequences",
    ):
        batch = unique_sequences[start:start + batch_size]
        outputs = batch_embed(model, batch, config)

        for sequence, output in zip(batch, outputs):
            sequence_embeddings[sequence] = _pool_embedding(output)

    # Return one embedding row per observation, even when sequences repeat.
    X = np.stack([sequence_embeddings[seq] for seq in sequences])
    y = pd.to_numeric(data[target_col]).to_numpy(dtype=float)

    embeddings_df = pd.DataFrame(
        X,
        columns=[f"emb_{i}" for i in range(X.shape[1])],
    )
    embeddings_df["label"] = y
    embeddings_df["sequence"] = sequences

    metadata_cols = metadata_cols or []
    for col in metadata_cols:
        if col in data.columns:
            embeddings_df[col] = data[col].to_numpy()

    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        embeddings_df.to_csv(output_path, index=False)

    return X, y, embeddings_df

def load_embeddings(file_path: str):
    """Load cached embeddings without treating metadata as model features."""
    df = pd.read_csv(file_path)

    feature_cols = [
        col for col in df.columns
        if col.startswith("emb_") and col[4:].isdigit()
    ]

    if not feature_cols:
        raise ValueError(f"No embedding columns found in {file_path!r}")

    if "label" not in df.columns or "sequence" not in df.columns:
        raise ValueError("Embedding file must contain 'label' and 'sequence'.")

    X = df[feature_cols].to_numpy(dtype=np.float32)
    y = pd.to_numeric(df["label"]).to_numpy(dtype=float)

    return X, y, df

"""Biological interaction graph construction.

Constructs confidence-weighted biological graphs G_B = (V, E, W) and Laplacians L_B = D - W
from sources such as STRING (protein-protein interactions) and pathway databases.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# STRING API base URL (v12)
_STRING_API = "https://string-db.org/api/json/network"

# Default STRING combined score threshold (400 = medium confidence)
DEFAULT_STRING_THRESHOLD = 400


# ---------------------------------------------------------------------------
# STRING-based graph
# ---------------------------------------------------------------------------

def fetch_string_interactions(
    gene_list: list[str],
    species: int = 9606,
    threshold: int = DEFAULT_STRING_THRESHOLD,
    cache_path: Optional[str | Path] = None,
) -> pd.DataFrame:
    """Fetch protein-protein interactions from STRING API.

    Parameters
    ----------
    gene_list : list[str]
        List of gene symbols (e.g., ["TP53", "BRCA1", "EGFR"]).
    species : int
        NCBI taxonomy ID. 9606 = human, 10090 = mouse.
    threshold : int
        Minimum combined interaction score (0–1000). 400 = medium confidence.
    cache_path : str or Path or None
        If provided, cache the API response to this path (avoids repeat calls).

    Returns
    -------
    pandas.DataFrame
        Columns: ``preferredName_A``, ``preferredName_B``, ``score``.
    """
    import pandas as pd
    import requests

    if cache_path is not None:
        cache_path = Path(cache_path)
        if cache_path.exists():
            logger.info("Loading STRING interactions from cache: %s", cache_path)
            return pd.read_csv(cache_path)

    logger.info(
        "Fetching STRING interactions for %d genes (species=%d, threshold=%d)",
        len(gene_list), species, threshold,
    )

    # STRING API: POST with gene identifiers
    params = {
        "identifiers": "\r".join(gene_list),
        "species": species,
        "required_score": threshold,
        "caller_identity": "mopa_framework",
    }

    resp = requests.post(_STRING_API, data=params, timeout=120)
    resp.raise_for_status()

    data = resp.json()
    if not data:
        logger.warning("STRING returned 0 interactions for the provided gene list.")
        return pd.DataFrame(columns=["preferredName_A", "preferredName_B", "score"])

    df = pd.DataFrame(data)
    logger.info("  STRING returned %d interactions", len(df))

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache_path, index=False)
        logger.info("  Cached to %s", cache_path)

    return df[["preferredName_A", "preferredName_B", "score"]]


def build_string_graph(
    gene_list: list[str],
    species: int = 9606,
    threshold: int = DEFAULT_STRING_THRESHOLD,
    cache_path: Optional[str | Path] = None,
) -> tuple[sp.csr_matrix, np.ndarray, list[str]]:
    """Build a confidence-weighted adjacency matrix from STRING interactions.

    Implements W = A ⊙ C (Hadamard product of adjacency × confidence matrix).
    STRING combined score ∈ [0, 1000] is normalized to C ∈ [0, 1].

    Parameters
    ----------
    gene_list : list[str]
        Gene symbols to build graph over.
    species : int
        NCBI taxonomy ID.
    threshold : int
        Minimum STRING confidence score.
    cache_path : str or Path or None
        Cache path for the STRING response CSV.

    Returns
    -------
    W : scipy.sparse.csr_matrix
        Confidence-weighted symmetric adjacency matrix (n × n).
    C : np.ndarray
        Confidence matrix (dense, n × n), values in [0, 1].
    gene_order : list[str]
        Ordered list of genes corresponding to matrix rows/columns.
    """

    df = fetch_string_interactions(gene_list, species, threshold, cache_path)

    gene_set = set(gene_list)
    gene_order = sorted(gene_set)
    gene_to_idx = {g: i for i, g in enumerate(gene_order)}
    n = len(gene_order)

    logger.info("Building STRING graph: %d genes", n)

    rows, cols, data = [], [], []
    for _, row in df.iterrows():
        g_a = row["preferredName_A"]
        g_b = row["preferredName_B"]
        score = float(row["score"]) / 1000.0  # normalize to [0, 1]

        if g_a in gene_to_idx and g_b in gene_to_idx:
            i, j = gene_to_idx[g_a], gene_to_idx[g_b]
            # Symmetric
            rows.extend([i, j])
            cols.extend([j, i])
            data.extend([score, score])

    W = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    C = W.toarray()

    logger.info("  Graph: %d nodes, %d non-zero edges", n, W.nnz // 2)
    return W, C, gene_order


# ---------------------------------------------------------------------------
# Graph Laplacian
# ---------------------------------------------------------------------------

def build_graph_laplacian(
    W: sp.csr_matrix,
    normalize: bool = False,
) -> sp.csr_matrix:
    """Compute the graph Laplacian L_B = D - W.

    The Laplacian is positive semi-definite for non-negative symmetric W.

    Parameters
    ----------
    W : scipy.sparse.csr_matrix
        Confidence-weighted symmetric adjacency matrix.
    normalize : bool
        If True, return the symmetric normalized Laplacian:
        L_sym = D^{-1/2} L D^{-1/2} = I - D^{-1/2} W D^{-1/2}

    Returns
    -------
    scipy.sparse.csr_matrix
        Graph Laplacian L_B ≽ 0.
    """
    n = W.shape[0]
    degree = np.asarray(W.sum(axis=1)).flatten()
    D = sp.diags(degree, format="csr")
    L = D - W

    if normalize:
        # Symmetric normalized: L_sym = I - D^{-1/2} W D^{-1/2}
        d_inv_sqrt = np.where(degree > 0, 1.0 / np.sqrt(degree), 0.0)
        D_inv_sqrt = sp.diags(d_inv_sqrt, format="csr")
        L = sp.eye(n) - D_inv_sqrt @ W @ D_inv_sqrt

    logger.debug("Laplacian: shape=%s, nnz=%d", L.shape, L.nnz)
    return L.tocsr()


# ---------------------------------------------------------------------------
# Save / Load
# ---------------------------------------------------------------------------

def save_graph(
    W: sp.csr_matrix,
    gene_order: list[str],
    path: str | Path,
) -> None:
    """Save the biological graph to a compressed NPZ archive.

    Parameters
    ----------
    W : scipy.sparse.csr_matrix
        Confidence-weighted adjacency.
    gene_order : list[str]
        Ordered gene names.
    path : str or Path
        Output path (will be created with .npz extension if needed).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sp.save_npz(str(path), W)
    # Save gene order alongside
    gene_path = path.with_suffix(".genes.txt")
    gene_path.write_text("\n".join(gene_order), encoding="utf-8")

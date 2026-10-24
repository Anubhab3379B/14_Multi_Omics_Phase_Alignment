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

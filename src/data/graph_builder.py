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

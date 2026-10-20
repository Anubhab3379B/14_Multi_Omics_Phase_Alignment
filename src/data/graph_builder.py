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


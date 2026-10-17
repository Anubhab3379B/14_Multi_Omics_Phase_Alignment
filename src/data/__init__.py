"""
src/data/__init__.py — Public API for the MOPA data ingestion layer.

Exposes the core loader functions and MuData preprocessing routines
for all supported multi-omics datasets.
"""
from src.data.loader import (
    load_neurips_cite,
    load_neurips_multiome,
    load_10x_h5_multiome,
    load_mouse_gastrulation_multiome,
    load_tar_dataset,
)
from src.data.preprocess import (
    preprocess_rna,
    preprocess_atac,
    preprocess_protein,
    preprocess_methylation,
    build_mudata,
)

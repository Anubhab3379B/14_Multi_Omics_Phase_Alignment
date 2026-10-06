# Multi-Omics Phase Alignment (MOPA)

[![Tests](https://img.shields.io/badge/pytest-21%20passed-brightgreen.svg)](tests/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**MOPA** is a deep generative systems biology framework for continuous stochastic trajectory reconstruction and multi-clock phase alignment across single-cell multi-omics modalities (scRNA-seq, scATAC-seq, and CITE-seq).

---

## Key Features

1. **Modality-Specific Decoders**: Biophysically aligned emission distributions (Negative Binomial for overdispersed counts, Bernoulli for binary chromatin peaks, Gaussian for normalized protein levels).
2. **Precision-Weighted Multimodal Fusion**: Bayesian Product-of-Experts fusion that weights input channels dynamically by their inverse variance.
3. **Biological Graph Regularization**: Graph Laplacian Dirichlet energy constraints derived from molecular interaction networks (STRING, pathway topologies).
4. **Neural Latent SDE**: Continuous-time state evolution modeling deterministic drift (lineage commitment) and state-dependent diffusion (stochastic fate branching).
5. **Monotonic Neural CDF Time-Warping**: Learnable, strictly increasing functions $\tau_m(t)$ capturing asynchronous progression rates across modalities without temporal inversion.
6. **Entropic Optimal Transport**: Differentiable cross-modal distribution alignment via Sinkhorn divergence.

---

## Benchmark Results

Evaluated across **18 systematic ablation experiments** ($M_0$ through $M_5$) on three paired single-cell cohorts:

- **NeurIPS 2021 CITE-seq** (10,000 cells; RNA + Protein)
- **NeurIPS 2021 Multiome** (5,000 cells; RNA + ATAC)
- **10x Genomics PBMC Multiome** (10,629 cells; RNA + ATAC)

### Summary Performance Table

| Dataset | Model | Smoothness ($\mathcal{S}_L \downarrow$) | Spearman ($\rho_t \uparrow$) | FOSCTTM ($\downarrow$) | Total Loss ($\downarrow$) |
|:---|:---|:---:|:---:|:---:|:---:|
| **NeurIPS CITE-seq** | $M_0$ (Baseline VAE) | $2.2646 \times 10^{-4}$ | 0.1512 | 0.1203 | 0.6551 |
| NeurIPS CITE-seq | $M_1$ (+Graph) | $2.0390 \times 10^{-4}$ | 0.2947 | 0.1080 | 0.6566 |
| NeurIPS CITE-seq | $M_2$ (+SDE) | $1.8313 \times 10^{-4}$ | 0.0383 | 0.1103 | 0.6571 |
| NeurIPS CITE-seq | $M_3$ (+OT) | $1.8719 \times 10^{-4}$ | 0.1598 | 0.1167 | 0.6580 |
| NeurIPS CITE-seq | $M_4$ (+Phase) | $1.8705 \times 10^{-4}$ | 0.1383 | 0.5404 | 0.6573 |
| NeurIPS CITE-seq | **$M_5$ (Full MOPA)** | **$1.7504 \times 10^{-4}$ (-22.6%)** | **0.1837** | 0.5634 | 0.6577 |
| **NeurIPS Multiome** | $M_0$ (Baseline VAE) | $3.6822 \times 10^{-4}$ | 0.0083 | 0.6448 | 0.4024 |
| NeurIPS Multiome | $M_3$ (+OT) | $3.2055 \times 10^{-4}$ | 0.0005 | 0.5517 | 0.4025 |
| NeurIPS Multiome | **$M_5$ (Full MOPA)** | **$3.3049 \times 10^{-4}$** | **0.0068** | 0.5542 | 0.4041 |
| **10x PBMC Multiome**| $M_0$ (Baseline VAE) | $9.4988 \times 10^{-5}$ | 0.0001 | 0.5108 | 0.2680 |
| 10x PBMC Multiome| $M_3$ (+OT) | $8.4691 \times 10^{-5}$ | 0.0018 | 0.5123 | 0.2683 |
| 10x PBMC Multiome| **$M_5$ (Full MOPA)** | **$8.6016 \times 10^{-5}$** | **0.0122** | 0.5042 | 0.2678 |

---

## Installation & Setup

### Environment Setup

```bash
# Clone the repository
git clone https://github.com/your-org/14_Multi_Omics_Phase_Alignment.git
cd 14_Multi_Omics_Phase_Alignment

# Install dependencies
pip install -r requirements.txt
```

### Running Tests

```bash
pytest
```

---

## Reproducibility Workflow

To reproduce all results, metrics, figures, and manuscripts from scratch:

```bash
# 1. Run all 18 ablation experiments
python scripts/run_all_experiments.py

# 2. Compute comprehensive benchmark evaluation metrics
python scripts/evaluate_extended_metrics.py

# 3. Generate high-resolution publication figures (300 DPI)
python scripts/generate_publication_figures.py

# 4. Compile styled Word manuscript
python scripts/generate_docx_manuscript.py
```

---

## Repository Structure

```text
├── configs/               # Hyperparameter and model configuration files
├── data/
│   ├── raw/               # Downloaded source datasets
│   └── processed/         # Standardized MuData (.h5mu) files
├── docs/                  # Technical documentation and experiment ledgers
├── results/
│   ├── experiments/       # Serialized metrics and training histories
│   ├── figures/           # Rendered 300 DPI publication plots
│   └── checkpoints/       # Saved PyTorch model weights
├── scripts/               # Automation scripts for training and evaluation
├── specifications/        # Formal mathematical and architecture specs
├── src/
│   ├── data/              # Ingestion, preprocessing, graph building, splitting
│   ├── losses/            # ELBO, SDE, Sinkhorn OT, Graph Dirichlet losses
│   ├── models/            # Encoders, fusion, SDE, time-warp, phase detection
│   └── training/          # Training loops, schedulers, and evaluation metrics
├── tests/                 # Pytest test suite (unit and pipeline tests)
└── thesis/                # Manuscript markdown source
```

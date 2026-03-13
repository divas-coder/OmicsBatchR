# OmicsBatchR

> A unified R package for batch correction, benchmarking, and comparison across single-cell and bulk RNA-seq data

**Built by the Divya Mishra**

---

## Overview

Batch effects are one of the biggest challenges in omics data analysis. OmicsBatchR solves this by giving you **one function to run 7 batch correction methods**, **automatic benchmarking**, and **publication-ready comparison plots** — all in R.

```r
results <- run_batch_correction(seurat_obj, batch_key = "batch", methods = "all")
best    <- select_best_method(results, batch_key = "batch")
plot_batch_comparison(results, color_by = "batch")
```

---

## Supported methods

| Method | Language | Type | Best for |
|--------|----------|------|----------|
| **Harmony** | R | Graph-based | Large datasets, fast |
| **Seurat CCA** | R | Anchor-based | Multi-dataset integration |
| **ComBat-seq** | R | Statistical model | Count-level bulk RNA-seq |
| **FastMNN** | R | Mutual nearest neighbours | Single-cell, conservative |
| **Scanorama** | Python | Panoramic stitching | Multi-batch panoramic |
| **BBKNN** | Python | Graph-based | Very large atlases |
| **scVI** | Python | Deep learning (VAE) | Complex batch structures |

---

## Installation

```r
# Install from GitHub
devtools::install_github("panomiq/OmicsBatchR")

# For Python-based methods (scVI, Scanorama, BBKNN)
reticulate::conda_create("r-reticulate", python_version = "3.9")
reticulate::conda_install("r-reticulate",
  packages = c("scvi-tools", "scanorama", "bbknn", "scanpy", "anndata"),
  pip = TRUE)
```

---

## Quick start

```r
library(OmicsBatchR)

# Step 1: Run all batch correction methods
results <- run_batch_correction(
  seurat_obj = your_seurat,
  batch_key  = "batch",
  methods    = "all",
  n_dims     = 30
)

# Step 2: Benchmark methods with LISI and kBET scores
scores <- benchmark_methods(
  results,
  batch_key    = "batch",
  celltype_key = "cell_type"
)

# Step 3: Auto-select the best method
best <- select_best_method(results, batch_key = "batch")
cat("Best method:", best$best_method)

# Step 4: Visualise UMAP comparison across all methods
plot_batch_comparison(results, color_by = c("batch", "cell_type"))

# Step 5: Plot benchmark scores
plot_benchmark_scores(scores)
```

---

## Benchmarking metrics

| Metric | Meaning | Direction |
|--------|---------|-----------|
| **iLISI** | Integration LISI — measures batch mixing | Higher = better |
| **cLISI** | Cell-type LISI — measures biology preservation | Lower = better |
| **kBET acceptance rate** | Local neighbourhood batch test | Higher = better |
| **Overall score** | Weighted composite (40% iLISI + 30% cLISI + 30% kBET) | Higher = better |

---

## Output

Running `run_batch_correction()` returns a list of corrected Seurat objects, one per method:

```
results
├── harmony      → Seurat obj with "harmony" + "umap_harmony" reductions
├── seurat_cca   → Seurat obj with "pca_cca" + "umap_seurat_cca" reductions
├── combat       → Seurat obj with "pca_combat" + "umap_combat" reductions
├── fastmnn      → Seurat obj with "mnn" + "umap_fastmnn" reductions
├── scanorama    → Seurat obj with "scanorama" + "umap_scanorama" reductions
├── bbknn        → Seurat obj with "bbknn_pca" + "umap_bbknn" reductions
└── scvi         → Seurat obj with "scvi" + "umap_scvi" reductions
```

---

## Package structure

```
OmicsBatchR/
├── R/
│   ├── OmicsBatchR-package.R      # Package documentation + imports
│   ├── run_batch_correction.R     # Main function + all method implementations
│   ├── benchmark_methods.R        # LISI, kBET scoring + auto-selection
│   └── plot_comparisons.R         # UMAP grid + benchmark bar charts
├── tests/testthat/
│   └── test-batch-correction.R    # Unit tests
├── vignettes/
│   └── OmicsBatchR.Rmd            # Full tutorial
├── DESCRIPTION
└── README.md
```

---



Contact:  [Upwork](https://www.upwork.com/freelancers/~01a3c803aa6c9d4418) | [LinkedIn](https://www.linkedin.com/in/dr-divya-mishra/)

---

## Citation

```
OmicsBatchR: Comprehensive
Batch Correction for Single-Cell and Bulk Omics Data.
GitHub: https://github.com/OmicsBatchR
```

## License

MIT License

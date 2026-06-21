# Primitives — authoritative spec

This document is filled in by Wave 2-A. The frozen interface lives in `packages/primitives/lattice_primitives/_interface.py`. Every registered primitive must appear here with its inputs, outputs, params, citation key, and sanity checks.

See AGENTS.md "How to add a primitive" for the workflow.

## Categories

- `qc` — quality control metrics, doublet detection
- `preprocess` — filtering, normalization, HVG selection
- `dim_reduce` — PCA, UMAP, t-SNE
- `cluster` — Leiden, Louvain
- `integrate` — Harmony, scVI, BBKNN
- `annotate` — CellTypist, mLLMCelltype consensus
- `de` — differential expression, pseudobulk + DESeq2
- `plot` — figures with publication defaults

## Registered primitives

### calculate_qc_metrics

| Field | Value |
|---|---|
| **Name** | `calculate_qc_metrics` |
| **Category** | `qc` |
| **Scanpy version** | `>=1.11.0,<1.12` |
| **Citation key** | `scanpy_calculate_qc_metrics` |
| **Row-modifying** | No |
| **File** | `lattice_primitives/qc/calculate_qc_metrics.py` |

**Requires:**
- No prerequisite obs columns
- No required layer_state (runs on raw data)

**Produces (obs columns):**
- `total_counts` — total UMI count per cell
- `n_genes_by_counts` — number of genes detected per cell
- `pct_counts_mt` — percentage of counts from mitochondrial genes

**Produces (var columns):**
- `mean_counts` — mean UMI count across cells
- `n_cells_by_counts` — number of cells in which the gene is detected

**Immutable obs columns:** `condition` (case/control swap protection)

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `mt_prefix` | `str` | `"MT-"` | Prefix for mitochondrial gene names. `"MT-"` for human; `"mt-"` for mouse. Auto-detects lowercase variant if no genes match. |

**Sanity checks:**
- `warn` if median `pct_counts_mt` > 20% — possibly low-quality tissue or stressed cells

---

### filter_cells_min_counts

| Field | Value |
|---|---|
| **Name** | `filter_cells_min_counts` |
| **Category** | `preprocess` |
| **Scanpy version** | `>=1.11.0,<1.12` |
| **Citation key** | `scanpy_filter_cells` |
| **Row-modifying** | Yes (removes cells) |
| **File** | `lattice_primitives/preprocess/filter_cells_basic.py` |

**Requires:**
- `total_counts` in obs (from `calculate_qc_metrics`)

**Produces:**
- Cells with `total_counts < min_counts_per_cell` are removed
- `condition` values of surviving cells are preserved bit-for-bit

**Immutable obs columns:** `condition` (checked per-cell, subset semantics for row-modifying primitives)

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `min_counts_per_cell` | `int` | `500` | Minimum total UMI count for a cell to be retained. Cells below this threshold are removed as low-quality. |

**Sanity checks:**
- `warn` if >50% of cells removed
- `error` if >80% of cells removed (threshold almost certainly too aggressive)

---

### filter_genes_min_cells

| Field | Value |
|---|---|
| **Name** | `filter_genes_min_cells` |
| **Category** | `preprocess` |
| **Scanpy version** | `>=1.11.0,<1.12` |
| **Citation key** | `scanpy_filter_cells` |
| **Row-modifying** | No (removes genes/var, not cells/obs) |
| **File** | `lattice_primitives/preprocess/filter_cells_basic.py` |

**Requires:**
- No required obs columns

**Produces:**
- Genes detected in fewer than `min_cells_per_gene` cells are removed from `adata.var`
- `n_vars` decreases; `n_obs` is unchanged

**Immutable obs columns:** `condition`

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `min_cells_per_gene` | `int` | `3` | Minimum number of cells a gene must be detected in to be retained. |

**Sanity checks:**
- `warn` if >80% of genes removed (unusually high gene dropout)

---

### normalize_total_log1p

| Field | Value |
|---|---|
| **Name** | `normalize_total_log1p` |
| **Category** | `preprocess` |
| **Scanpy version** | `>=1.11.0,<1.12` |
| **Citation key** | `theis_normalization_2019` |
| **Row-modifying** | No |
| **File** | `lattice_primitives/preprocess/normalize_total_log1p.py` |

**Requires:**
- `layer_state = "raw_counts"` (set by the ingest primitive)

**Produces:**
- `layer_state = "log_normalized"`
- `adata.layers["counts_normalized"]` — pre-normalization raw count copy
- `adata.X` replaced with `log1p(normalized counts)`

**Immutable obs columns:** `condition` (case/control swap protection)

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `target_sum` | `float \| None` | `None` | Per-cell normalization target. `None` uses the median total count across cells (Theis lab recommendation). `1e4` = CP10k, `1e6` = CPM. |
| `exclude_highly_expressed` | `bool` | `False` | Exclude top-1% most expressed genes from the normalization denominator. Useful for erythrocyte-rich or highly imbalanced tissues. |

**Sanity checks:**
- `error` if `max(adata.X) > 15` after log1p — input may not have been raw integer counts
- `warn` if median normalized total per cell deviates from `log1p(target_sum)` by > 1.0

---

## Citation keys

All citation keys reference entries in `lattice_primitives/citations.yaml`.

| Key | Reference |
|---|---|
| `scanpy_calculate_qc_metrics` | Wolf et al., SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 2018. DOI: 10.1186/s13059-017-1382-0 |
| `scanpy_filter_cells` | Wolf et al., SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 2018. DOI: 10.1186/s13059-017-1382-0 |
| `theis_normalization_2019` | Luecken & Theis, Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology* 2019. DOI: 10.15252/msb.20188746 |

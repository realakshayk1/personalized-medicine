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

_Populated in Wave 2-A._

"""Tests for the neighbors and umap primitives."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
import lattice_primitives.dim_reduce.umap  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_pca(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    return out


def test_neighbors_builds_graph(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    adata_out, _, _ = run_primitive("neighbors", adata, n_neighbors=15, random_state=0)
    assert "neighbors" in adata_out.uns
    assert "connectivities" in adata_out.obsp


def test_neighbors_requires_pca(adata_log_normalized: AnnData) -> None:
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    # No PCA run -> in-function guard raises ValueError
    with pytest.raises(ValueError, match="PCA embedding"):
        run_primitive("neighbors", out, n_neighbors=15)


def test_umap_adds_embedding(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    adata, _, _ = run_primitive("neighbors", adata, n_neighbors=15, random_state=0)
    adata_out, _, _ = run_primitive("umap", adata, random_state=0)
    assert "X_umap" in adata_out.obsm
    assert adata_out.obsm["X_umap"].shape == (adata_out.n_obs, 2)
    assert np.all(np.isfinite(adata_out.obsm["X_umap"]))


def test_umap_requires_neighbors(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    with pytest.raises(ValueError, match="neighbor graph"):
        run_primitive("umap", adata, random_state=0)

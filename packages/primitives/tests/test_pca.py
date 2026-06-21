"""Tests for the pca primitive."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
from lattice_primitives import run_primitive
from lattice_primitives.registry import SchemaViolation


def _hvg(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    return out


def test_happy_path_adds_embedding(adata_log_normalized: AnnData) -> None:
    adata = _hvg(adata_log_normalized)
    adata_out, _, result = run_primitive("pca", adata, n_comps=30, random_state=0)
    assert "X_pca" in adata_out.obsm
    assert adata_out.obsm["X_pca"].shape == (adata_out.n_obs, 30)
    assert "pca" in adata_out.uns
    assert np.all(np.isfinite(adata_out.obsm["X_pca"]))
    assert result.primitive == "pca"


def test_rejects_without_hvg(adata_log_normalized: AnnData) -> None:
    # pca requires the highly_variable var column
    with pytest.raises(SchemaViolation):
        run_primitive("pca", adata_log_normalized, n_comps=30)


def test_reproducible_with_seed(adata_log_normalized: AnnData) -> None:
    a1 = _hvg(adata_log_normalized.copy())
    a2 = _hvg(adata_log_normalized.copy())
    out1, _, _ = run_primitive("pca", a1, n_comps=20, random_state=7)
    out2, _, _ = run_primitive("pca", a2, n_comps=20, random_state=7)
    np.testing.assert_allclose(
        np.abs(out1.obsm["X_pca"]), np.abs(out2.obsm["X_pca"]), rtol=1e-5, atol=1e-5
    )

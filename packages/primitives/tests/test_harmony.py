"""Tests for the harmony batch-integration primitive."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
import lattice_primitives.integrate.harmony  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_pca(adata: AnnData) -> AnnData:
    """Run highly_variable_genes + pca via real run_primitive."""
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    return out


def _add_two_batches(adata: AnnData) -> AnnData:
    """Assign two alternating batch labels in adata.obs['batch']."""
    labels = ["b1" if i % 2 == 0 else "b2" for i in range(adata.n_obs)]
    adata.obs["batch"] = labels
    return adata


def _add_single_batch(adata: AnnData) -> AnnData:
    """Assign a single batch label (integration is a no-op)."""
    adata.obs["batch"] = ["b1"] * adata.n_obs
    return adata


def test_harmony_adds_embedding(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    adata = _add_two_batches(adata)
    adata_out, warnings, _ = run_primitive("harmony", adata, batch_key="batch", random_state=0)
    assert "X_pca_harmony" in adata_out.obsm
    x_harmony = np.asarray(adata_out.obsm["X_pca_harmony"])
    assert x_harmony.shape == np.asarray(adata_out.obsm["X_pca"]).shape
    assert np.all(np.isfinite(x_harmony))
    # Two real batches: no single-batch no-op warning expected.
    assert not any("one batch level" in w.message for w in warnings)


def test_harmony_requires_pca(adata_log_normalized: AnnData) -> None:
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    out = _add_two_batches(out)
    # No PCA run -> in-function guard raises ValueError.
    with pytest.raises(ValueError, match="PCA embedding"):
        run_primitive("harmony", out, batch_key="batch", random_state=0)


def test_harmony_requires_batch_key(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    # batch_key column absent -> in-function guard raises ValueError.
    with pytest.raises(ValueError, match="batch_key"):
        run_primitive("harmony", adata, batch_key="nonexistent", random_state=0)


def test_harmony_warns_on_single_batch(adata_log_normalized: AnnData) -> None:
    adata = _through_pca(adata_log_normalized)
    adata = _add_single_batch(adata)
    _, warnings, _ = run_primitive("harmony", adata, batch_key="batch", random_state=0)
    assert any(w.severity == "warn" and "one batch level" in w.message for w in warnings)


def test_harmony_reproducible(adata_log_normalized: AnnData) -> None:
    # Two independent copies from the same clean input so each run is identical.
    adata1 = _add_two_batches(_through_pca(adata_log_normalized.copy()))
    out1, _, _ = run_primitive("harmony", adata1, batch_key="batch", random_state=0)
    emb1 = np.asarray(out1.obsm["X_pca_harmony"]).copy()

    adata2 = _add_two_batches(_through_pca(adata_log_normalized.copy()))
    out2, _, _ = run_primitive("harmony", adata2, batch_key="batch", random_state=0)
    emb2 = np.asarray(out2.obsm["X_pca_harmony"])

    assert np.allclose(emb1, emb2)

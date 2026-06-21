"""Tests for the leiden clustering primitive."""

from __future__ import annotations

import pytest
from anndata import AnnData

import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_neighbors(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    return out


def test_leiden_finds_two_populations(adata_log_normalized: AnnData) -> None:
    adata = _through_neighbors(adata_log_normalized)
    adata_out, _, result = run_primitive("leiden", adata, resolution=1.0, random_state=0)
    assert "leiden" in adata_out.obs.columns
    n_clusters = adata_out.obs["leiden"].nunique()
    # Two separable populations in the fixture; expect >= 2 clusters.
    assert n_clusters >= 2
    assert result.primitive == "leiden"


def test_leiden_requires_neighbors(adata_log_normalized: AnnData) -> None:
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    with pytest.raises(ValueError, match="neighbor graph"):
        run_primitive("leiden", out, resolution=1.0)


def test_leiden_reproducible(adata_log_normalized: AnnData) -> None:
    a1 = _through_neighbors(adata_log_normalized.copy())
    a2 = _through_neighbors(adata_log_normalized.copy())
    out1, _, _ = run_primitive("leiden", a1, resolution=1.0, random_state=42)
    out2, _, _ = run_primitive("leiden", a2, resolution=1.0, random_state=42)
    assert list(out1.obs["leiden"]) == list(out2.obs["leiden"])


def test_leiden_preserves_condition(adata_log_normalized: AnnData) -> None:
    adata = _through_neighbors(adata_log_normalized)
    before = list(adata.obs["condition"])
    out, _, _ = run_primitive("leiden", adata, resolution=1.0, random_state=0)
    # immutable_obs_columns guard: condition must be untouched
    assert list(out.obs["condition"]) == before

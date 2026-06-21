"""Tests for the rank_genes_groups (marker DE) primitive."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.de.rank_genes_groups  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_leiden(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    return out


def test_rank_genes_groups_happy_path(adata_log_normalized: AnnData) -> None:
    adata = _through_leiden(adata_log_normalized)
    adata_out, warnings, result = run_primitive(
        "rank_genes_groups", adata, groupby="leiden", method="wilcoxon"
    )
    assert "rank_genes_groups" in adata_out.uns
    rgg = adata_out.uns["rank_genes_groups"]
    # Scanpy structured result has per-group columns.
    assert len(rgg["names"].dtype.names) >= 2
    assert result.primitive == "rank_genes_groups"
    # No error-severity warnings on clean data.
    assert all(w.severity != "error" for w in warnings)


def test_rank_genes_groups_rejects_missing_groupby(adata_log_normalized: AnnData) -> None:
    adata = _through_leiden(adata_log_normalized)
    with pytest.raises(ValueError, match="groupby column"):
        run_primitive("rank_genes_groups", adata, groupby="does_not_exist")


def test_rank_genes_groups_preserves_condition(adata_log_normalized: AnnData) -> None:
    adata = _through_leiden(adata_log_normalized)
    before = list(adata.obs["condition"])
    out, _, _ = run_primitive("rank_genes_groups", adata, groupby="leiden")
    assert list(out.obs["condition"]) == before


def test_rank_genes_groups_rejects_singleton_group(adata_log_normalized: AnnData) -> None:
    """A group of size 1 is rejected with a clear error.

    scanpy itself raises on singleton groups (one-vs-rest is statistically
    undefined), so the primitive guards before the call and surfaces an
    actionable message rather than scanpy's opaque one.
    """
    adata = _through_leiden(adata_log_normalized)
    # Inject a singleton cluster into a fresh categorical column.
    labels = np.array(["A"] * adata.n_obs, dtype=object)
    labels[0] = "singleton"
    adata.obs["custom_groups"] = labels
    adata.obs["custom_groups"] = adata.obs["custom_groups"].astype("category")
    with pytest.raises(ValueError, match="fewer than 2 cells"):
        run_primitive("rank_genes_groups", adata, groupby="custom_groups", method="wilcoxon")

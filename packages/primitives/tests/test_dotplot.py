"""Tests for the dotplot primitive (figure emission)."""

from __future__ import annotations

import pytest
from anndata import AnnData

import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.de.rank_genes_groups  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
import lattice_primitives.plot.dotplot  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_markers(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    out, _, _ = run_primitive("rank_genes_groups", out, groupby="leiden", method="wilcoxon")
    return out


def test_dotplot_emits_figure(adata_log_normalized: AnnData) -> None:
    adata = _through_markers(adata_log_normalized)
    adata_out, _, result = run_primitive("dotplot", adata, groupby="leiden", n_genes=5)
    assert len(result.figures) == 1
    assert result.figures[0].startswith("data:image/png;base64,")
    assert "_lattice_figures" not in adata_out.uns


def test_dotplot_requires_rank_genes_groups(adata_log_normalized: AnnData) -> None:
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    with pytest.raises(ValueError, match="rank_genes_groups"):
        run_primitive("dotplot", out, groupby="leiden")


def test_dotplot_rejects_missing_groupby(adata_log_normalized: AnnData) -> None:
    adata = _through_markers(adata_log_normalized)
    with pytest.raises(ValueError, match="not found in adata.obs"):
        run_primitive("dotplot", adata, groupby="does_not_exist")

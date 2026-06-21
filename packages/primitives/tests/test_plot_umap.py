"""Tests for the plot_umap primitive (figure emission)."""

from __future__ import annotations

import pytest
from anndata import AnnData

import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
import lattice_primitives.dim_reduce.umap  # noqa: F401  (register)
import lattice_primitives.plot.plot_umap  # noqa: F401  (register)
from lattice_primitives import run_primitive


def _through_leiden(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("umap", out, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    return out


def test_plot_umap_emits_figure(adata_log_normalized: AnnData) -> None:
    adata = _through_leiden(adata_log_normalized)
    adata_out, _, result = run_primitive("plot_umap", adata, color="leiden")
    # Figure is lifted into the result and removed from .uns
    assert len(result.figures) == 1
    assert result.figures[0].startswith("data:image/png;base64,")
    assert "_lattice_figures" not in adata_out.uns


def test_plot_umap_requires_umap(adata_log_normalized: AnnData) -> None:
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    with pytest.raises(ValueError, match="UMAP embedding"):
        run_primitive("plot_umap", out, color="leiden")


def test_plot_umap_rejects_missing_color(adata_log_normalized: AnnData) -> None:
    adata = _through_leiden(adata_log_normalized)
    with pytest.raises(ValueError, match="not found in adata.obs"):
        run_primitive("plot_umap", adata, color="does_not_exist")

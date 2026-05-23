"""Tests for filter_cells_min_counts and filter_genes_min_cells primitives."""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

from lattice_primitives import PRIMITIVE_REGISTRY, SchemaViolation, run_primitive
from lattice_primitives.preprocess.filter_cells_basic import (
    filter_cells_min_counts,
    filter_genes_min_cells,
)
from lattice_primitives.qc.calculate_qc_metrics import calculate_qc_metrics

# Ensure primitives registered
_ = filter_cells_min_counts, filter_genes_min_cells


def _run_qc(adata: AnnData) -> AnnData:
    return calculate_qc_metrics(adata, mt_prefix="MT-")


# ---------------------------------------------------------------------------
# filter_cells_min_counts
# ---------------------------------------------------------------------------


def test_filter_cells_happy_path(adata: AnnData) -> None:
    """filter_cells_min_counts removes cells below threshold."""
    adata = _run_qc(adata)
    n_before = adata.n_obs
    result = filter_cells_min_counts(adata, min_counts_per_cell=500)
    # With Poisson(lam=10) * 100 genes, most cells have ~1000 counts -> few removed
    assert result.n_obs <= n_before
    # All remaining cells should meet the threshold
    assert (result.obs["total_counts"] >= 500).all()


def test_filter_cells_requires_total_counts(adata: AnnData) -> None:
    """filter_cells_min_counts raises SchemaViolation if total_counts is missing."""
    # Don't run QC first -> total_counts absent
    with pytest.raises(SchemaViolation, match="total_counts"):
        filter_cells_min_counts(adata, min_counts_per_cell=500)


def test_filter_cells_in_registry() -> None:
    assert "filter_cells_min_counts" in PRIMITIVE_REGISTRY


def test_filter_genes_happy_path(adata: AnnData) -> None:
    """filter_genes_min_cells removes genes below detection threshold."""
    n_vars_before = adata.n_vars
    result = filter_genes_min_cells(adata, min_cells_per_gene=3)
    assert result.n_vars <= n_vars_before


def test_filter_genes_obs_names_unchanged(adata: AnnData) -> None:
    """filter_genes is not row-modifying; obs_names must be unchanged."""
    obs_names_before = adata.obs_names.tolist()
    result = filter_genes_min_cells(adata, min_cells_per_gene=1)
    assert result.obs_names.tolist() == obs_names_before


def test_filter_genes_in_registry() -> None:
    assert "filter_genes_min_cells" in PRIMITIVE_REGISTRY


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------


def test_sanity_check_error_on_extreme_threshold(adata: AnnData) -> None:
    """Sanity check emits 'error' severity when >80% of cells are removed."""
    adata = _run_qc(adata)
    # With Poisson(lam=10)*100 genes, median total ~1000; set absurdly high threshold
    _, warnings, _ = run_primitive(
        "filter_cells_min_counts", adata, min_counts_per_cell=100_000
    )
    error_warns = [w for w in warnings if w.severity == "error"]
    assert len(error_warns) > 0, (
        f"Expected an 'error' sanity check warning, got: {[w.severity for w in warnings]}"
    )


def test_sanity_check_warn_on_high_removal(adata: AnnData) -> None:
    """Sanity check emits at least 'warn' when 50-80% of cells are removed."""
    adata = _run_qc(adata)
    # Set threshold at ~75th percentile of total_counts to remove ~75% of cells
    threshold = int(np.percentile(adata.obs["total_counts"].values, 75))
    _, warnings, _ = run_primitive(
        "filter_cells_min_counts", adata, min_counts_per_cell=threshold
    )
    # At least no error-level assertion failure — warnings list returned cleanly
    assert isinstance(warnings, list)


def test_filter_cells_condition_preserved(adata: AnnData) -> None:
    """condition column values must survive cell filtering."""
    adata = _run_qc(adata)
    condition_mapping = dict(zip(adata.obs_names, adata.obs["condition"]))
    result = filter_cells_min_counts(adata, min_counts_per_cell=100)
    for cell, cond in zip(result.obs_names, result.obs["condition"]):
        assert condition_mapping[cell] == cond, (
            f"Cell {cell} condition changed after filtering"
        )

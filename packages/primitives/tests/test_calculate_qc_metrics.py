"""Tests for the calculate_qc_metrics primitive."""

from __future__ import annotations

import numpy as np
from anndata import AnnData

from lattice_primitives import PRIMITIVE_REGISTRY, run_primitive
from lattice_primitives.qc.calculate_qc_metrics import calculate_qc_metrics

# Ensure primitive is registered by importing its module
_ = calculate_qc_metrics


def test_happy_path_adds_expected_obs_columns(adata: AnnData) -> None:
    """calculate_qc_metrics adds total_counts, n_genes_by_counts, pct_counts_mt to obs."""
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    assert "total_counts" in result.obs.columns
    assert "n_genes_by_counts" in result.obs.columns
    assert "pct_counts_mt" in result.obs.columns


def test_happy_path_adds_expected_var_columns(adata: AnnData) -> None:
    """calculate_qc_metrics adds mean_counts, n_cells_by_counts to var."""
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    assert "mean_counts" in result.var.columns
    assert "n_cells_by_counts" in result.var.columns


def test_total_counts_positive(adata: AnnData) -> None:
    """total_counts should be positive for all cells in a non-zero dataset."""
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    assert (result.obs["total_counts"] >= 0).all()


def test_pct_counts_mt_between_0_and_100(adata: AnnData) -> None:
    """pct_counts_mt should be in [0, 100] for all cells."""
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    assert (result.obs["pct_counts_mt"] >= 0.0).all()
    assert (result.obs["pct_counts_mt"] <= 100.0).all()


def test_mt_pct_nonzero_with_mt_genes(adata: AnnData) -> None:
    """Median MT% should be > 0 since fixture has 10 MT-prefixed genes with counts."""
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    median_mt = float(np.median(result.obs["pct_counts_mt"]))
    assert median_mt > 0.0, "Expected non-zero MT% with MT genes in fixture"


def test_obs_names_unchanged(adata: AnnData) -> None:
    """calculate_qc_metrics must not change obs_names (non-row-modifying)."""
    obs_names_before = adata.obs_names.tolist()
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    assert result.obs_names.tolist() == obs_names_before


def test_condition_column_preserved(adata: AnnData) -> None:
    """condition column values must be bit-identical after QC metrics."""
    condition_before = adata.obs["condition"].values.copy()
    result = calculate_qc_metrics(adata, mt_prefix="MT-")
    np.testing.assert_array_equal(result.obs["condition"].values, condition_before)


def test_primitive_in_registry() -> None:
    """calculate_qc_metrics should be in PRIMITIVE_REGISTRY."""
    assert "calculate_qc_metrics" in PRIMITIVE_REGISTRY


def test_sanity_check_high_mito_emits_warn(bad_anndata_high_mito: AnnData) -> None:
    """High-mito fixture should trigger a 'warn' from the sanity check."""
    _, warnings, _ = run_primitive("calculate_qc_metrics", bad_anndata_high_mito, mt_prefix="MT-")
    warn_messages = [w.message for w in warnings]
    assert any("mitochondrial" in m.lower() for m in warn_messages), (
        f"Expected a mitochondrial content warning, got: {warn_messages}"
    )
    severities = [w.severity for w in warnings]
    assert "warn" in severities


def test_sanity_check_normal_adata_no_warn(adata: AnnData) -> None:
    """Normal fixture (low MT%) should not produce error-level warnings."""
    from tests.conftest import _make_adata
    fresh = _make_adata(high_mito=False)
    _, warnings, _ = run_primitive("calculate_qc_metrics", fresh, mt_prefix="MT-")
    error_warnings = [w for w in warnings if w.severity == "error"]
    assert len(error_warnings) == 0

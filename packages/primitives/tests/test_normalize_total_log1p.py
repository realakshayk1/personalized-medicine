"""Tests for normalize_total_log1p primitive."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from lattice_primitives import PRIMITIVE_REGISTRY, SchemaViolation, run_primitive
from lattice_primitives.preprocess.normalize_total_log1p import normalize_total_log1p

# Ensure primitive is registered
_ = normalize_total_log1p


def test_normalize_happy_path_layer_state(adata_raw_counts_state: AnnData) -> None:
    """normalize_total_log1p sets layer_state to log_normalized."""
    result = normalize_total_log1p(adata_raw_counts_state)
    assert result.uns.get("_lattice_layer_state") == "log_normalized"


def test_normalize_adds_counts_normalized_layer(adata_raw_counts_state: AnnData) -> None:
    """counts_normalized layer must be added."""
    result = normalize_total_log1p(adata_raw_counts_state)
    assert "counts_normalized" in result.layers


def test_normalize_x_is_log_transformed(adata_raw_counts_state: AnnData) -> None:
    """After log1p, max value should be reasonable (< 15 for raw counts)."""
    result = normalize_total_log1p(adata_raw_counts_state)
    X = result.X
    if scipy.sparse.issparse(X):
        max_val = float(X.max())
    else:
        max_val = float(np.asarray(X).max())
    assert max_val < 15.0, f"Max value {max_val} after log1p suggests non-raw input"


def test_normalize_obs_names_unchanged(adata_raw_counts_state: AnnData) -> None:
    """normalize is not row-modifying; obs_names must be unchanged."""
    obs_before = adata_raw_counts_state.obs_names.tolist()
    result = normalize_total_log1p(adata_raw_counts_state)
    assert result.obs_names.tolist() == obs_before


def test_normalize_in_registry() -> None:
    assert "normalize_total_log1p" in PRIMITIVE_REGISTRY


def test_normalize_requires_raw_counts_state(adata: AnnData) -> None:
    """Should raise SchemaViolation if layer_state != raw_counts."""
    adata.uns["_lattice_layer_state"] = "log_normalized"  # wrong state
    with pytest.raises(SchemaViolation, match="layer_state"):
        normalize_total_log1p(adata)


def test_normalize_without_any_layer_state_raises(adata: AnnData) -> None:
    """Without _lattice_layer_state in uns, requires.layer_state='raw_counts' should fail."""
    assert "_lattice_layer_state" not in adata.uns
    with pytest.raises(SchemaViolation, match="layer_state"):
        normalize_total_log1p(adata)


def test_normalize_condition_preserved(adata_raw_counts_state: AnnData) -> None:
    """condition column must be bit-identical after normalization."""
    condition_before = adata_raw_counts_state.obs["condition"].values.copy()
    result = normalize_total_log1p(adata_raw_counts_state)
    np.testing.assert_array_equal(result.obs["condition"].values, condition_before)


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


def test_sanity_check_error_on_pre_logged_data() -> None:
    """Sanity check emits 'error' when input looks already log-transformed."""
    # Create adata with very large raw values so that after normalization
    # some cells/genes will have extremely high normalized counts -> log1p > 15
    n_obs, n_vars = 50, 20
    counts = np.ones((n_obs, n_vars), dtype=np.float32) * 0.1
    # One gene with massive count in one cell -> after normalize to median, it remains huge
    counts[:, 0] = 1e8  # massive outlier gene
    X = scipy.sparse.csr_matrix(counts)
    adata = AnnData(X=X)
    adata.obs_names = [f"C{i}" for i in range(n_obs)]
    adata.var_names = [f"G{i}" for i in range(n_vars)]
    adata.obs["condition"] = ["ctrl"] * n_obs
    adata.uns["_lattice_layer_state"] = "raw_counts"

    _, warnings, _ = run_primitive("normalize_total_log1p", adata)
    error_warns = [w for w in warnings if w.severity == "error"]
    assert len(error_warns) > 0, (
        f"Expected sanity check error for extreme values, got: {[w.message for w in warnings]}"
    )

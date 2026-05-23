"""Tests for the primitive registry mechanics."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

from lattice_primitives import (
    PRIMITIVE_REGISTRY,
    AnnDataState,
    ParamSpec,
    PrimitiveNotFound,
    SchemaViolation,
    primitive,
    run_primitive,
)

# ---------------------------------------------------------------------------
# Dummy primitive for registry tests
# ---------------------------------------------------------------------------


@primitive(
    name="test_dummy_noop",
    category="qc",
    requires=AnnDataState(obs_columns=["test_col"]),
    produces=AnnDataState(obs_columns=["result_col"]),
    params={
        "multiplier": ParamSpec(**{"type": "float", "default": 1.0, "description": "test param"}),
    },
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_calculate_qc_metrics",
    row_modifying=False,
)
def _test_dummy_noop(adata: AnnData, multiplier: float = 1.0) -> AnnData:
    adata.obs["result_col"] = adata.obs["test_col"] * multiplier
    return adata


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_registry_contains_dummy() -> None:
    """@primitive decorator registers the primitive by name."""
    assert "test_dummy_noop" in PRIMITIVE_REGISTRY


def test_registry_entry_fields() -> None:
    """Registry entry has correct metadata."""
    entry = PRIMITIVE_REGISTRY["test_dummy_noop"]
    assert entry.name == "test_dummy_noop"
    assert entry.category == "qc"
    assert entry.citation_key == "scanpy_calculate_qc_metrics"
    assert "multiplier" in entry.params


def test_dummy_has_primitive_name_attr() -> None:
    """Wrapped function has __primitive_name__ attribute."""
    assert _test_dummy_noop.__primitive_name__ == "test_dummy_noop"  # type: ignore[attr-defined]


def test_dummy_happy_path(adata: AnnData) -> None:
    """Primitive runs successfully when requires are met."""
    adata.obs["test_col"] = np.arange(adata.n_obs, dtype=float)
    result = _test_dummy_noop(adata, multiplier=2.0)
    assert "result_col" in result.obs.columns
    expected = np.arange(adata.n_obs, dtype=float) * 2.0
    np.testing.assert_array_almost_equal(result.obs["result_col"].values, expected)


def test_requires_missing_column_raises(adata: AnnData) -> None:
    """SchemaViolation is raised when required obs column is absent."""
    # 'test_col' not in adata.obs -> should raise
    if "test_col" in adata.obs.columns:
        del adata.obs["test_col"]
    with pytest.raises(SchemaViolation, match="test_col"):
        _test_dummy_noop(adata)


def test_run_primitive_not_found() -> None:
    """run_primitive raises PrimitiveNotFound for unknown name."""
    X = scipy.sparse.csr_matrix(np.zeros((5, 5)))
    adata = AnnData(X=X)
    with pytest.raises(PrimitiveNotFound):
        run_primitive("this_primitive_does_not_exist", adata)


def test_run_primitive_returns_result_structure(adata: AnnData) -> None:
    """run_primitive returns (AnnData, list[Warning_], PrimitiveResult)."""
    from lattice_primitives import PrimitiveResult
    adata.obs["test_col"] = np.arange(adata.n_obs, dtype=float)
    out_adata, warnings, result = run_primitive("test_dummy_noop", adata, multiplier=1.5)
    assert isinstance(result, PrimitiveResult)
    assert result.primitive == "test_dummy_noop"
    assert result.input_hash != ""
    assert result.output_hash != ""
    assert result.duration_sec >= 0.0
    assert isinstance(warnings, list)

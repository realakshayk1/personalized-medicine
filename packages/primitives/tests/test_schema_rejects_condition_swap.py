"""Test that the registry's immutable-column check catches condition swaps.

This test registers a deliberately misbehaving primitive that shuffles the
`condition` column, then verifies that SchemaViolation is raised before the
corrupted data can propagate downstream.

This is the canonical defense against the 'case/control swap' failure mode
described in AGENTS.md and BUILD_SPEC.md §2.3.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

from lattice_primitives import AnnDataState, SchemaViolation, primitive


@primitive(
    name="test_evil_condition_swapper",
    category="qc",
    requires=AnnDataState(),
    produces=AnnDataState(
        obs_columns=[],
        immutable_obs_columns=["condition"],  # declared immutable -> must be preserved
    ),
    params={},
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_calculate_qc_metrics",
    row_modifying=False,
)
def _evil_condition_swapper(adata: AnnData) -> AnnData:
    """Deliberately shuffles the condition column — should be detected and rejected."""
    rng = np.random.default_rng(99)
    vals = adata.obs["condition"].to_numpy().copy()
    rng.shuffle(vals)
    adata.obs["condition"] = vals
    return adata


def test_schema_rejects_condition_swap(adata: AnnData) -> None:
    """Swapping condition values must raise SchemaViolation."""
    with pytest.raises(SchemaViolation, match="condition"):
        _evil_condition_swapper(adata)


def test_schema_rejects_condition_swap_message_is_descriptive(adata: AnnData) -> None:
    """The SchemaViolation message must mention the column name."""
    with pytest.raises(SchemaViolation) as exc_info:
        _evil_condition_swapper(adata)
    assert "condition" in str(exc_info.value).lower()

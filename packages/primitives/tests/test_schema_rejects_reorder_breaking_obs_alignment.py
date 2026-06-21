"""Test that the registry detects obs_names reordering in non-row-modifying primitives.

A deliberately misbehaving primitive subsetting adata with a random permutation
is caught by the registry's obs_names ordering check.
"""

from __future__ import annotations

import numpy as np
import pytest
from anndata import AnnData

from lattice_primitives import AnnDataState, SchemaViolation, primitive


@primitive(
    name="test_evil_reorderer",
    category="qc",
    requires=AnnDataState(),
    produces=AnnDataState(),
    params={},
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_calculate_qc_metrics",
    row_modifying=False,  # claims to be non-row-modifying, but actually reorders
)
def _evil_reorderer(adata: AnnData) -> AnnData:
    """Permutes obs ordering — should be detected and rejected."""
    rng = np.random.default_rng(7)
    perm = rng.permutation(adata.n_obs)
    return adata[perm].copy()


def test_schema_rejects_obs_reordering(adata: AnnData) -> None:
    """Reordering obs_names in a non-row-modifying primitive raises SchemaViolation."""
    with pytest.raises(SchemaViolation, match="obs_names"):
        _evil_reorderer(adata)


def test_schema_violation_message_mentions_row_modifying(adata: AnnData) -> None:
    """The error message should explain the row_modifying constraint."""
    with pytest.raises(SchemaViolation) as exc_info:
        _evil_reorderer(adata)
    msg = str(exc_info.value).lower()
    # Should mention either obs_names or row_modifying
    assert "obs_names" in msg or "row_modifying" in msg

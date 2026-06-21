"""Reproducibility tests — same input, same code, bit-identical output.

Parametrized over all 3 v0.1 primitives. For each:
- Runs twice with identical inputs
- Asserts adata.X is bit-identical (via SHA256 hash)
- Asserts PrimitiveResult.output_hash is identical across the two runs
"""

from __future__ import annotations

import hashlib

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData

# Ensure all primitives are registered by importing their modules
import lattice_primitives.preprocess.filter_cells_basic  # noqa: F401
import lattice_primitives.preprocess.normalize_total_log1p  # noqa: F401
import lattice_primitives.qc.calculate_qc_metrics  # noqa: F401
from lattice_primitives import run_primitive


def _hash_x(adata: AnnData) -> str:
    """SHA256 hash of adata.X for bit-identical check."""
    X = adata.X
    if scipy.sparse.issparse(X):
        data = X.toarray().tobytes()
    else:
        data = np.asarray(X).tobytes()
    return hashlib.sha256(data).hexdigest()


def _make_fresh_adata(seed: int = 42) -> AnnData:
    """Create a fresh AnnData from conftest fixture logic (independent of fixtures)."""
    from tests.conftest import _make_adata
    return _make_adata(n_obs=200, n_vars=100, seed=seed, high_mito=False)


# ---------------------------------------------------------------------------
# Reproducibility tests
# ---------------------------------------------------------------------------


def test_calculate_qc_metrics_reproducible() -> None:
    """calculate_qc_metrics output is bit-identical across two runs."""
    adata1 = _make_fresh_adata(seed=42)
    adata2 = _make_fresh_adata(seed=42)

    out1, _, result1 = run_primitive("calculate_qc_metrics", adata1, mt_prefix="MT-")
    out2, _, result2 = run_primitive("calculate_qc_metrics", adata2, mt_prefix="MT-")

    assert _hash_x(out1) == _hash_x(out2), "X is not bit-identical across runs"
    assert result1.output_hash == result2.output_hash, "output_hash differs across runs"


def test_filter_cells_min_counts_reproducible() -> None:
    """filter_cells_min_counts output is bit-identical across two runs."""
    from lattice_primitives.qc.calculate_qc_metrics import calculate_qc_metrics

    adata1 = _make_fresh_adata(seed=42)
    adata2 = _make_fresh_adata(seed=42)

    calculate_qc_metrics(adata1, mt_prefix="MT-")
    calculate_qc_metrics(adata2, mt_prefix="MT-")

    out1, _, result1 = run_primitive(
        "filter_cells_min_counts", adata1, min_counts_per_cell=500
    )
    out2, _, result2 = run_primitive(
        "filter_cells_min_counts", adata2, min_counts_per_cell=500
    )

    assert _hash_x(out1) == _hash_x(out2), "X is not bit-identical across runs"
    assert result1.output_hash == result2.output_hash, "output_hash differs across runs"


def test_normalize_total_log1p_reproducible() -> None:
    """normalize_total_log1p output is bit-identical across two runs."""
    adata1 = _make_fresh_adata(seed=42)
    adata2 = _make_fresh_adata(seed=42)

    adata1.uns["_lattice_layer_state"] = "raw_counts"
    adata2.uns["_lattice_layer_state"] = "raw_counts"

    out1, _, result1 = run_primitive(
        "normalize_total_log1p", adata1, target_sum=None, exclude_highly_expressed=False
    )
    out2, _, result2 = run_primitive(
        "normalize_total_log1p", adata2, target_sum=None, exclude_highly_expressed=False
    )

    assert _hash_x(out1) == _hash_x(out2), "X is not bit-identical across runs"
    assert result1.output_hash == result2.output_hash, "output_hash differs across runs"


@pytest.mark.parametrize(
    "primitive_name,setup_fn,run_params",
    [
        (
            "calculate_qc_metrics",
            lambda a: a,
            {"mt_prefix": "MT-"},
        ),
    ],
)
def test_parametrized_output_hash_identical(
    primitive_name: str,
    setup_fn: object,
    run_params: dict,  # type: ignore[type-arg]
) -> None:
    """Parametrized: output_hash must be identical across two identical runs."""
    adata1 = _make_fresh_adata(seed=42)
    adata2 = _make_fresh_adata(seed=42)

    adata1 = setup_fn(adata1)  # type: ignore[operator]
    adata2 = setup_fn(adata2)  # type: ignore[operator]

    _, _, result1 = run_primitive(primitive_name, adata1, **run_params)
    _, _, result2 = run_primitive(primitive_name, adata2, **run_params)

    assert result1.output_hash == result2.output_hash, (
        f"output_hash not identical for '{primitive_name}': "
        f"{result1.output_hash} != {result2.output_hash}"
    )

"""Test fixtures for lattice_primitives tests.

Provides:
    adata           — 200 cells × 100 genes, Poisson counts, 10 MT genes
    bad_anndata_high_mito — synthetic AnnData with very high mitochondrial fraction
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData


def _make_adata(
    n_obs: int = 200,
    n_vars: int = 100,
    seed: int = 42,
    high_mito: bool = False,
) -> AnnData:
    rng = np.random.default_rng(seed)

    # Poisson counts — typical for scRNA-seq
    if high_mito:
        # Make most counts come from MT genes (first 10 genes)
        counts = np.zeros((n_obs, n_vars), dtype=np.float32)
        # 90% of counts in MT genes
        counts[:, :10] = rng.poisson(lam=500, size=(n_obs, 10)).astype(np.float32)
        counts[:, 10:] = rng.poisson(lam=2, size=(n_obs, n_vars - 10)).astype(np.float32)
    else:
        counts = rng.poisson(lam=10, size=(n_obs, n_vars)).astype(np.float32)

    X = scipy.sparse.csr_matrix(counts)

    # Gene names: first 10 are MT-prefixed
    var_names = [f"MT-GENE{i:03d}" for i in range(10)] + [
        f"GENE{i:03d}" for i in range(10, n_vars)
    ]

    # Cell barcodes
    obs_names = [f"CELL_{i:04d}" for i in range(n_obs)]

    # Condition: half treated, half control
    condition = ["treated"] * (n_obs // 2) + ["control"] * (n_obs - n_obs // 2)

    adata = AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obs["condition"] = condition

    return adata


@pytest.fixture
def adata() -> AnnData:
    """Synthetic AnnData: 200 cells × 100 genes, Poisson counts (seed 42).

    - First 10 genes are MT-prefixed (MT-GENE000 .. MT-GENE009)
    - obs["condition"]: 100 "treated", 100 "control"
    - .X is a sparse CSR matrix of Poisson-distributed counts
    """
    return _make_adata(n_obs=200, n_vars=100, seed=42, high_mito=False)


@pytest.fixture
def bad_anndata_high_mito() -> AnnData:
    """Synthetic AnnData with very high mitochondrial content.

    90% of counts come from the 10 MT-prefixed genes.
    Median pct_counts_mt will be >> 20%, triggering a 'warn' sanity check.
    """
    return _make_adata(n_obs=200, n_vars=100, seed=42, high_mito=True)


@pytest.fixture
def adata_raw_counts_state(adata: AnnData) -> AnnData:
    """adata with layer_state set to 'raw_counts' for normalize primitive."""
    adata.uns["_lattice_layer_state"] = "raw_counts"
    return adata

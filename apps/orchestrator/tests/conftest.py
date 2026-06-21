"""Shared test fixtures for the orchestrator test suite."""

from __future__ import annotations

import json
from collections.abc import AsyncIterator
from typing import Any

import anndata
import lattice_primitives.preprocess.filter_cells_basic  # noqa: F401
import lattice_primitives.preprocess.normalize_total_log1p  # noqa: F401

# ---------------------------------------------------------------------------
# Primitive registration — ensure modules are imported before any test runs
# ---------------------------------------------------------------------------
import lattice_primitives.qc.calculate_qc_metrics  # noqa: F401
import numpy as np
import pytest
import pytest_asyncio
import scipy.sparse
from lattice.models import Plan, PrimitiveInvocation
from lattice.provenance import ProvenanceStore
from lattice.sandbox import FakeSandbox

# ---------------------------------------------------------------------------
# Synthetic AnnData fixture (mirrors the primitives test fixture)
# ---------------------------------------------------------------------------


def make_synthetic_adata(
    n_obs: int = 50,
    n_vars: int = 100,
    with_condition: bool = True,
) -> anndata.AnnData:
    """Create a synthetic AnnData with raw integer counts."""
    import pandas as pd
    rng = np.random.default_rng(42)
    X = scipy.sparse.csr_matrix(
        rng.negative_binomial(5, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    )
    adata = anndata.AnnData(X=X)
    adata.obs_names = pd.Index([f"cell_{i}" for i in range(n_obs)])
    adata.var_names = pd.Index([f"gene_{i}" for i in range(n_vars)])
    if with_condition:
        adata.obs["condition"] = ["treated"] * (n_obs // 2) + ["control"] * (n_obs - n_obs // 2)
    # Mark layer_state as raw_counts so normalize_total_log1p is happy
    adata.uns["_lattice_layer_state"] = "raw_counts"
    return adata


@pytest.fixture
def synthetic_adata() -> anndata.AnnData:
    return make_synthetic_adata()


# ---------------------------------------------------------------------------
# FakeSandbox fixture
# ---------------------------------------------------------------------------


@pytest_asyncio.fixture
async def fake_sandbox(synthetic_adata: anndata.AnnData) -> AsyncIterator[FakeSandbox]:
    sb = FakeSandbox()
    await sb.upload_anndata(synthetic_adata)
    yield sb
    await sb.close()


# ---------------------------------------------------------------------------
# ProvenanceStore fixture (in-memory SQLite)
# ---------------------------------------------------------------------------


@pytest_asyncio.fixture
async def prov_store(tmp_path: Any) -> AsyncIterator[ProvenanceStore]:
    db_path = str(tmp_path / "test_prov.db")
    store = ProvenanceStore(db_path=db_path)
    await store.open()
    yield store
    await store.close()


# ---------------------------------------------------------------------------
# Stub chat function for Planner tests
# ---------------------------------------------------------------------------

_STUB_PLAN_4_STEPS: dict[str, Any] = {
    "session_id": "test-session",
    "steps": [
        {
            "primitive": "calculate_qc_metrics",
            "params": {"mt_prefix": "MT-"},
            "rationale": "Compute QC metrics.",
        },
        {
            "primitive": "filter_cells_min_counts",
            "params": {"min_counts_per_cell": 200},
            "rationale": "Filter low-count cells.",
        },
        {
            "primitive": "filter_genes_min_cells",
            "params": {"min_cells_per_gene": 3},
            "rationale": "Filter low-expression genes.",
        },
        {
            "primitive": "normalize_total_log1p",
            "params": {"target_sum": None, "exclude_highly_expressed": False},
            "rationale": "Normalize counts.",
        },
    ],
    "rationale": "Standard 4-step scRNA-seq pipeline.",
    "workflow": "standard_scrnaseq",
}


async def stub_chat_fn_4steps(system_prompt: str, user_message: str) -> str:
    """Stub that always returns the 4-step standard plan JSON."""
    return json.dumps(_STUB_PLAN_4_STEPS)


async def stub_chat_fn_unknown_primitive(system_prompt: str, user_message: str) -> str:
    """Stub that returns a plan with an unregistered primitive."""
    bad_plan = {**_STUB_PLAN_4_STEPS}
    bad_plan["steps"] = [
        {
            "primitive": "hallucinated_scanpy_func",
            "params": {},
            "rationale": "This primitive does not exist.",
        }
    ]
    return json.dumps(bad_plan)


@pytest.fixture
def stub_chat_4steps() -> Any:
    return stub_chat_fn_4steps


@pytest.fixture
def stub_chat_unknown() -> Any:
    return stub_chat_fn_unknown_primitive


# ---------------------------------------------------------------------------
# 4-step plan fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def standard_4step_plan() -> Plan:
    return Plan(
        session_id="test-session",
        steps=[
            PrimitiveInvocation(
                primitive="calculate_qc_metrics",
                params={"mt_prefix": "MT-"},
            ),
            PrimitiveInvocation(
                primitive="filter_cells_min_counts",
                params={"min_counts_per_cell": 200},
            ),
            PrimitiveInvocation(
                primitive="filter_genes_min_cells",
                params={"min_cells_per_gene": 3},
            ),
            PrimitiveInvocation(
                primitive="normalize_total_log1p",
                params={"target_sum": None, "exclude_highly_expressed": False},
            ),
        ],
        rationale="Standard 4-step scRNA-seq pipeline.",
        workflow="standard_scrnaseq",
    )


# ---------------------------------------------------------------------------
# Minimal h5ad bytes fixture for upload tests
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_h5ad_bytes() -> bytes:
    """Create a minimal h5ad file in memory and return its bytes."""
    import tempfile
    anndata.settings.allow_write_nullable_strings = True
    adata = make_synthetic_adata(n_obs=10, n_vars=20, with_condition=False)
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
        tmp_path = f.name
    adata.write_h5ad(tmp_path)
    with open(tmp_path, "rb") as f:
        return f.read()

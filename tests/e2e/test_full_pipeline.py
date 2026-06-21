"""End-to-end smoke: synthetic AnnData → Plan → Executor → Provenance.

Drives the orchestrator's full execution path with a stubbed planner LLM and
FakeSandbox (in-process primitives). No network calls. Validates that:

  - the planner can be bypassed with a hand-crafted Plan that passes validation,
  - the executor walks the plan, yields the expected ExecuteEvent sequence,
  - every step lands in the provenance store with input/output hashes,
  - the immutable_obs check on `condition` survives the full chain.
"""

from __future__ import annotations

import anndata

# Import primitive modules so @primitive decorators register them.
import lattice_primitives.preprocess.filter_cells_basic  # noqa: F401
import lattice_primitives.preprocess.normalize_total_log1p  # noqa: F401
import lattice_primitives.qc.calculate_qc_metrics  # noqa: F401
import numpy as np
import pytest
import scipy.sparse
from lattice.executor import run_plan
from lattice.models import (
    Plan,
    PlanCompletedEvent,
    PrimitiveInvocation,
    StepCompletedEvent,
    StepStartedEvent,
)
from lattice.provenance import ProvenanceStore
from lattice.sandbox import FakeSandbox


def _synthetic_anndata(seed: int = 42) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    n_obs, n_vars = 200, 100
    counts = rng.poisson(0.5, size=(n_obs, n_vars)).astype(np.float32)
    X = scipy.sparse.csr_matrix(counts)
    var_names = [f"MT-{i}" if i < 10 else f"GENE{i}" for i in range(n_vars)]
    adata = anndata.AnnData(
        X=X,
        obs={"condition": ["treated"] * (n_obs // 2) + ["control"] * (n_obs // 2)},
        var={"gene_symbols": var_names},
    )
    adata.var_names = var_names
    adata.uns["_lattice_layer_state"] = "raw_counts"
    return adata


@pytest.mark.asyncio
async def test_full_pipeline_end_to_end(tmp_path):
    db_path = tmp_path / "provenance.sqlite"
    store = ProvenanceStore(str(db_path))
    await store.open()
    try:
        session_id = "e2e-test-session"
        await store.create_session(
            session_id=session_id,
            started_at="2026-05-23T00:00:00Z",
            filename="synthetic.h5ad",
            sha256="x" * 64,
            n_obs=200,
            n_vars=100,
        )

        sandbox = FakeSandbox()
        await sandbox.upload_anndata(_synthetic_anndata())

        plan = Plan(
            session_id=session_id,
            workflow="standard_scrnaseq",
            rationale="e2e test plan",
            steps=[
                PrimitiveInvocation(primitive="calculate_qc_metrics", params={}),
                PrimitiveInvocation(
                    primitive="filter_cells_min_counts",
                    params={"min_counts_per_cell": 5},
                ),
                PrimitiveInvocation(
                    primitive="filter_genes_min_cells",
                    params={"min_cells_per_gene": 3},
                ),
                PrimitiveInvocation(primitive="normalize_total_log1p", params={}),
            ],
        )

        events = [e async for e in run_plan(plan, sandbox, store)]

        # Event shape: 4 × (started, completed) + 1 × plan_completed
        starts = [e for e in events if isinstance(e, StepStartedEvent)]
        completes = [e for e in events if isinstance(e, StepCompletedEvent)]
        plan_done = [e for e in events if isinstance(e, PlanCompletedEvent)]
        assert len(starts) == 4
        assert len(completes) == 4
        assert len(plan_done) == 1
        assert plan_done[0].n_steps == 4

        # Provenance has 4 rows
        log = await store.get_provenance(session_id)
        assert len(log.steps) == 4
        assert [s.primitive for s in log.steps] == [
            "calculate_qc_metrics",
            "filter_cells_min_counts",
            "filter_genes_min_cells",
            "normalize_total_log1p",
        ]
        for step in log.steps:
            assert step.input_hash and step.output_hash
            assert step.duration_sec >= 0

        # Final AnnData state: log_normalized, condition column preserved
        final = await sandbox.download_anndata()
        assert final.uns["_lattice_layer_state"] == "log_normalized"
        assert "treated" in set(final.obs["condition"])
        assert "control" in set(final.obs["condition"])
        assert float(final.X.max()) < 15.0  # post-log1p sanity
    finally:
        await store.close()

"""Tests for the Executor (run_plan)."""

from __future__ import annotations

from typing import Any

import pytest
from lattice.executor import run_plan
from lattice.models import (
    Plan,
    PlanCompletedEvent,
    StepCompletedEvent,
    StepFailedEvent,
    StepStartedEvent,
)
from lattice.provenance import ProvenanceStore
from lattice.sandbox import FakeSandbox


class TestExecutor:
    @pytest.mark.asyncio
    async def test_4_events_in_order(
        self,
        standard_4step_plan: Plan,
        fake_sandbox: FakeSandbox,
        prov_store: ProvenanceStore,
    ) -> None:
        """4-step plan emits 4×(started+completed) + plan_completed = 9 events."""
        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        events = []
        async for event in run_plan(
            plan=standard_4step_plan,
            sandbox=fake_sandbox,
            store=prov_store,
        ):
            events.append(event)

        # Classify event types
        started = [e for e in events if isinstance(e, StepStartedEvent)]
        completed = [e for e in events if isinstance(e, StepCompletedEvent)]
        failed = [e for e in events if isinstance(e, StepFailedEvent)]
        plan_done = [e for e in events if isinstance(e, PlanCompletedEvent)]

        assert len(started) == 4
        assert len(completed) == 4
        assert len(failed) == 0
        assert len(plan_done) == 1
        assert plan_done[0].n_steps == 4

    @pytest.mark.asyncio
    async def test_events_in_correct_order(
        self,
        standard_4step_plan: Plan,
        fake_sandbox: FakeSandbox,
        prov_store: ProvenanceStore,
    ) -> None:
        """For each step: started then completed; plan_completed last."""
        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        events = []
        async for event in run_plan(
            plan=standard_4step_plan,
            sandbox=fake_sandbox,
            store=prov_store,
        ):
            events.append(event)

        # Verify order: started[0], completed[0], started[1], completed[1], ..., plan_completed
        expected_types = []
        for i in range(4):
            expected_types.append(("step_started", i))
            expected_types.append(("step_completed", i))
        expected_types.append(("plan_completed", None))

        actual_types = []
        for e in events:
            idx = getattr(e, "step_index", None)
            actual_types.append((e.type, idx))  # type: ignore[attr-defined]

        assert actual_types == expected_types

    @pytest.mark.asyncio
    async def test_provenance_has_4_rows(
        self,
        standard_4step_plan: Plan,
        fake_sandbox: FakeSandbox,
        prov_store: ProvenanceStore,
    ) -> None:
        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        async for _ in run_plan(
            plan=standard_4step_plan,
            sandbox=fake_sandbox,
            store=prov_store,
        ):
            pass

        prov = await prov_store.get_provenance("test-session")
        assert prov is not None
        assert len(prov.steps) == 4

    @pytest.mark.asyncio
    async def test_provenance_records_primitives(
        self,
        standard_4step_plan: Plan,
        fake_sandbox: FakeSandbox,
        prov_store: ProvenanceStore,
    ) -> None:
        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        async for _ in run_plan(
            plan=standard_4step_plan,
            sandbox=fake_sandbox,
            store=prov_store,
        ):
            pass

        prov = await prov_store.get_provenance("test-session")
        assert prov is not None
        primitive_names = [s.primitive for s in prov.steps]
        assert primitive_names == [
            "calculate_qc_metrics",
            "filter_cells_min_counts",
            "filter_genes_min_cells",
            "normalize_total_log1p",
        ]

    @pytest.mark.asyncio
    async def test_step_failed_halts_execution(
        self,
        prov_store: ProvenanceStore,
    ) -> None:
        """A failing step emits step_failed and stops; no further steps run."""


        # Create a plan with a bad second step (wrong param type)
        bad_plan = Plan(
            session_id="test-session",
            steps=[
                # First step succeeds
                __import__("lattice.models", fromlist=["PrimitiveInvocation"]).PrimitiveInvocation(
                    primitive="calculate_qc_metrics",
                    params={"mt_prefix": "MT-"},
                ),
                # Second step: intentionally call nonexistent primitive
                __import__("lattice.models", fromlist=["PrimitiveInvocation"]).PrimitiveInvocation(
                    primitive="calculate_qc_metrics",  # real but we'll patch sandbox
                    params={"mt_prefix": "MT-"},
                ),
            ],
            rationale="Test failure halting.",
            workflow="standard_scrnaseq",
        )

        # Use a sandbox that raises on second call
        call_count = 0

        class FailAfterFirstSandbox:
            async def upload_anndata(self, adata: Any) -> None:  # type: ignore[override]
                pass

            async def run_primitive(self, name: str, params: Any) -> Any:  # type: ignore[override]
                nonlocal call_count
                call_count += 1
                if call_count >= 2:
                    raise RuntimeError("Simulated sandbox failure")
                # Return a minimal PrimitiveResult for first call
                from lattice.models import PrimitiveResult
                return PrimitiveResult(
                    primitive=name,
                    params=params,
                    input_hash="abc",
                    output_hash="def",
                    duration_sec=0.01,
                )

            async def download_anndata(self) -> Any:  # type: ignore[override]
                raise NotImplementedError

            async def close(self) -> None:
                pass

        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        events = []
        async for event in run_plan(
            plan=bad_plan,
            sandbox=FailAfterFirstSandbox(),  # type: ignore[arg-type]
            store=prov_store,
        ):
            events.append(event)

        event_types = [e.type for e in events]  # type: ignore[attr-defined]
        assert "step_failed" in event_types
        # plan_completed should NOT be in events (halted)
        assert "plan_completed" not in event_types
        # Only 2 started events (step 0 started, step 1 started before fail)
        started = [e for e in events if e.type == "step_started"]  # type: ignore[attr-defined]
        assert len(started) == 2

    @pytest.mark.asyncio
    async def test_immutable_obs_survives(
        self,
        standard_4step_plan: Plan,
        fake_sandbox: FakeSandbox,
        prov_store: ProvenanceStore,
    ) -> None:
        """The 'condition' obs column is preserved bit-for-bit through all steps."""
        # Get initial condition values
        adata_before = await fake_sandbox.download_anndata()
        condition_before = adata_before.obs.get("condition")

        await prov_store.create_session(
            session_id="test-session",
            started_at="2026-05-22T00:00:00Z",
        )

        # Execute all steps
        events = []
        async for event in run_plan(
            plan=standard_4step_plan,
            sandbox=fake_sandbox,
            store=prov_store,
        ):
            events.append(event)

        adata_after = await fake_sandbox.download_anndata()

        # If condition existed before, it must survive
        if condition_before is not None:
            assert "condition" in adata_after.obs.columns
            # Values of surviving cells must match
            surviving_cells = set(adata_after.obs_names)
            orig_conditions = dict(zip(adata_before.obs_names, condition_before))
            for cell in surviving_cells:
                assert adata_after.obs.loc[cell, "condition"] == orig_conditions[cell]

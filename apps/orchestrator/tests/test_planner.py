"""Tests for the Planner class."""

from __future__ import annotations

from typing import Any

import pytest
from lattice.planner import Planner, UnknownPrimitiveError


class TestPlannerWithStubLLM:
    @pytest.mark.asyncio
    async def test_4step_plan_validates(self, stub_chat_4steps: Any) -> None:
        planner = Planner(chat_fn=stub_chat_4steps, workflow_name="standard_scrnaseq")
        plan = await planner.plan(
            session_id="test-session",
            user_message="Run standard QC and preprocessing",
            anndata_summary={"n_obs": 100, "n_vars": 200, "obs_columns": [], "filename": "test.h5ad"},
        )
        assert plan.session_id == "test-session"
        assert len(plan.steps) == 4
        primitives = [s.primitive for s in plan.steps]
        assert primitives == [
            "calculate_qc_metrics",
            "filter_cells_min_counts",
            "filter_genes_min_cells",
            "normalize_total_log1p",
        ]
        assert plan.workflow == "standard_scrnaseq"

    @pytest.mark.asyncio
    async def test_unknown_primitive_rejected(self, stub_chat_unknown: Any) -> None:
        planner = Planner(chat_fn=stub_chat_unknown, workflow_name="standard_scrnaseq")
        with pytest.raises(UnknownPrimitiveError, match="hallucinated_scanpy_func"):
            await planner.plan(
                session_id="test-session",
                user_message="Run something weird",
                anndata_summary={},
            )

    @pytest.mark.asyncio
    async def test_plan_params_preserved(self, stub_chat_4steps: Any) -> None:
        planner = Planner(chat_fn=stub_chat_4steps, workflow_name="standard_scrnaseq")
        plan = await planner.plan(
            session_id="sess-x",
            user_message="Standard pipeline",
            anndata_summary={"n_obs": 50, "n_vars": 100, "obs_columns": [], "filename": "x.h5ad"},
        )
        filter_step = plan.steps[1]
        assert filter_step.primitive == "filter_cells_min_counts"
        assert filter_step.params["min_counts_per_cell"] == 200

    @pytest.mark.asyncio
    async def test_session_id_injected(self, stub_chat_4steps: Any) -> None:
        """session_id is always overridden with the one passed to plan()."""
        planner = Planner(chat_fn=stub_chat_4steps, workflow_name="standard_scrnaseq")
        plan = await planner.plan(
            session_id="my-specific-session",
            user_message="Run QC",
            anndata_summary={},
        )
        assert plan.session_id == "my-specific-session"

    @pytest.mark.asyncio
    async def test_workflow_not_in_allowed_rejected(self) -> None:
        """A plan with a primitive outside workflow allowed_primitives is rejected."""
        import json
        bad_plan = {
            "session_id": "s",
            "steps": [
                {
                    "primitive": "leiden_clustering",  # not in standard_scrnaseq
                    "params": {},
                    "rationale": "Cluster cells.",
                }
            ],
            "rationale": "test",
            "workflow": "standard_scrnaseq",
        }

        async def bad_chat(system: str, user: str) -> str:
            return json.dumps(bad_plan)

        planner = Planner(chat_fn=bad_chat, workflow_name="standard_scrnaseq")
        with pytest.raises(UnknownPrimitiveError):
            await planner.plan(session_id="s", user_message="cluster", anndata_summary={})

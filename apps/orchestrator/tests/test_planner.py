"""Tests for the Planner class."""

from __future__ import annotations

import json
from typing import Any

import pytest
from lattice.models import Plan
from lattice.planner import (
    Planner,
    UnknownPrimitiveError,
    load_few_shot_examples,
    load_system_prompt,
    load_workflow,
    render_few_shot_block,
)


class TestPlannerWithStubLLM:
    @pytest.mark.asyncio
    async def test_4step_plan_validates(self, stub_chat_4steps: Any) -> None:
        planner = Planner(chat_fn=stub_chat_4steps, workflow_name="standard_scrnaseq")
        plan = await planner.plan(
            session_id="test-session",
            user_message="Run standard QC and preprocessing",
            anndata_summary={
                "n_obs": 100,
                "n_vars": 200,
                "obs_columns": [],
                "filename": "test.h5ad",
            },
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
            anndata_summary={
                "n_obs": 50,
                "n_vars": 100,
                "obs_columns": [],
                "filename": "x.h5ad",
            },
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
            await planner.plan(
                session_id="s", user_message="cluster", anndata_summary={}
            )

    @pytest.mark.asyncio
    async def test_valid_multistep_plan_passes(self) -> None:
        """A valid multi-step standard plan (QC -> ... -> leiden) validates."""
        steps = [
            {"primitive": "calculate_qc_metrics", "params": {"mt_prefix": "MT-"}},
            {
                "primitive": "filter_cells_min_counts",
                "params": {"min_counts_per_cell": 500},
            },
            {
                "primitive": "filter_genes_min_cells",
                "params": {"min_cells_per_gene": 3},
            },
            {"primitive": "normalize_total_log1p", "params": {"target_sum": None}},
            {"primitive": "highly_variable_genes", "params": {"n_top_genes": 2000}},
            {"primitive": "pca", "params": {"n_comps": 50, "random_state": 0}},
            {
                "primitive": "neighbors",
                "params": {"n_neighbors": 15, "random_state": 0},
            },
            {"primitive": "umap", "params": {"min_dist": 0.5, "random_state": 0}},
            {"primitive": "leiden", "params": {"resolution": 1.0, "random_state": 0}},
        ]
        good_plan = {
            "session_id": "s",
            "steps": steps,
            "rationale": "Standard QC + clustering.",
            "workflow": "standard_scrnaseq",
        }

        async def good_chat(system: str, user: str) -> str:
            return json.dumps(good_plan)

        planner = Planner(chat_fn=good_chat, workflow_name="standard_scrnaseq")
        plan = await planner.plan(
            session_id="my-sess",
            user_message="Run standard QC and clustering",
            anndata_summary={
                "n_obs": 2700,
                "n_vars": 32738,
                "obs_columns": [],
                "filename": "p.h5ad",
            },
        )
        assert len(plan.steps) == 9
        assert [s.primitive for s in plan.steps] == [s["primitive"] for s in steps]
        assert plan.session_id == "my-sess"

    @pytest.mark.asyncio
    async def test_cross_condition_de_plan_validates(self) -> None:
        """A pseudobulk -> pseudobulk_deseq2 plan validates under cross_condition_de."""
        de_plan = {
            "session_id": "s",
            "steps": [
                {
                    "primitive": "pseudobulk",
                    "params": {
                        "sample_col": "sample",
                        "groups_col": "cell_type",
                        "mode": "sum",
                    },
                },
                {
                    "primitive": "pseudobulk_deseq2",
                    "params": {
                        "condition_col": "condition",
                        "reference_level": "control",
                        "test_level": "treated",
                        "design": "~condition",
                    },
                },
            ],
            "rationale": "Cross-condition DE via pseudobulk + DESeq2.",
            "workflow": "cross_condition_de",
        }

        async def de_chat(system: str, user: str) -> str:
            return json.dumps(de_plan)

        planner = Planner(chat_fn=de_chat, workflow_name="cross_condition_de")
        plan = await planner.plan(
            session_id="de-sess",
            user_message="Compare treated vs control per cell type",
            anndata_summary={"obs_columns": ["sample", "condition", "cell_type"]},
        )
        assert [s.primitive for s in plan.steps] == ["pseudobulk", "pseudobulk_deseq2"]

    @pytest.mark.asyncio
    async def test_out_of_workflow_primitive_rejected(self) -> None:
        """A standard-pipeline primitive proposed under cross_condition_de is rejected."""
        bad_plan = {
            "session_id": "s",
            "steps": [{"primitive": "leiden", "params": {}, "rationale": "cluster"}],
            "rationale": "wrong workflow",
            "workflow": "cross_condition_de",
        }

        async def bad_chat(system: str, user: str) -> str:
            return json.dumps(bad_plan)

        planner = Planner(chat_fn=bad_chat, workflow_name="cross_condition_de")
        with pytest.raises(UnknownPrimitiveError):
            await planner.plan(
                session_id="s", user_message="cluster", anndata_summary={}
            )


# ---------------------------------------------------------------------------
# Few-shot examples + enriched system prompt
# ---------------------------------------------------------------------------


class TestFewShotExamples:
    def test_examples_file_loads_as_valid_jsonl(self) -> None:
        """planner_examples.jsonl loads and every line is a well-formed object."""
        examples = load_few_shot_examples()
        assert len(examples) >= 4, "expected at least 4 few-shot examples"
        for ex in examples:
            assert "user_message" in ex
            assert "plan" in ex
            assert "anndata_summary" in ex

    def test_every_example_plan_is_valid_against_its_workflow(self) -> None:
        """Each example Plan parses and uses only its workflow's allowed primitives."""
        from lattice_primitives import PRIMITIVE_REGISTRY  # noqa: F401
        import lattice_primitives.all_primitives  # noqa: F401

        allowed_by_workflow: dict[str, set[str]] = {}
        examples = load_few_shot_examples()
        assert examples
        for ex in examples:
            plan = Plan.model_validate(ex["plan"])
            wf = plan.workflow
            if wf not in allowed_by_workflow:
                allowed_by_workflow[wf] = set(
                    load_workflow(wf).get("allowed_primitives", [])
                )
            allowed = allowed_by_workflow[wf]
            assert plan.steps, f"example for '{wf}' has no steps"
            for step in plan.steps:
                assert step.primitive in PRIMITIVE_REGISTRY, (
                    f"'{step.primitive}' not registered (workflow {wf})"
                )
                assert step.primitive in allowed, (
                    f"'{step.primitive}' not in allowed_primitives for '{wf}'"
                )

    def test_examples_cover_required_scenarios(self) -> None:
        """The four required scenarios are represented in the few-shot set."""
        examples = load_few_shot_examples()
        plans = [Plan.model_validate(ex["plan"]) for ex in examples]
        prim_sets = [{s.primitive for s in p.steps} for p in plans]
        workflows = {p.workflow for p in plans}

        # (a) standard QC + clustering + UMAP
        assert any({"calculate_qc_metrics", "leiden", "umap"} <= ps for ps in prim_sets)
        # (b) full pipeline through markers + figures
        assert any(
            {"rank_genes_groups", "plot_umap", "dotplot", "heatmap", "volcano"} <= ps
            for ps in prim_sets
        )
        # (c) cross-condition DE
        assert "cross_condition_de" in workflows
        assert any({"pseudobulk", "pseudobulk_deseq2"} <= ps for ps in prim_sets)
        # (d) multi-batch case (per-batch HVG, integration mentioned)
        assert any("highly_variable_genes" in ps and "pca" in ps for ps in prim_sets)

    def test_render_block_empty_for_no_examples(self) -> None:
        assert render_few_shot_block([]) == ""

    def test_render_block_includes_example_content(self) -> None:
        examples = load_few_shot_examples()
        block = render_few_shot_block(examples)
        assert "Few-shot examples" in block
        assert "Plan:" in block
        # The block must contain at least one real primitive name.
        assert "calculate_qc_metrics" in block


class TestEnrichedSystemPrompt:
    def test_system_prompt_describes_primitives(self) -> None:
        """The system prompt now contains primitive descriptions + constraints."""
        prompt = load_system_prompt()
        for name in [
            "calculate_qc_metrics",
            "normalize_total_log1p",
            "highly_variable_genes",
            "pca",
            "neighbors",
            "umap",
            "leiden",
            "rank_genes_groups",
            "pseudobulk",
            "pseudobulk_deseq2",
        ]:
            assert name in prompt, f"system prompt missing description for {name}"
        assert "allowed_primitives" in prompt
        # ordering prerequisite language present
        assert "after" in prompt.lower()

    def test_planner_system_prompt_includes_examples(self) -> None:
        """Planner assembles system prompt + few-shot block (chat_fn seam intact)."""
        planner = Planner(workflow_name="standard_scrnaseq")
        assert "Few-shot examples" in planner._system_prompt
        assert "calculate_qc_metrics" in planner._system_prompt

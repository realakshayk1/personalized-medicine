"""Executor: walks a Plan and yields ExecuteEvents, recording provenance."""

from __future__ import annotations

import datetime
from collections.abc import AsyncIterator
from typing import Any

from loguru import logger

from lattice.models import (
    ExecuteEvent,
    Plan,
    PlanCompletedEvent,
    PrimitiveResult,
    StepCompletedEvent,
    StepFailedEvent,
    StepStartedEvent,
)
from lattice.provenance import ProvenanceStore
from lattice.sandbox import SandboxProtocol


async def run_plan(
    plan: Plan,
    sandbox: SandboxProtocol,
    store: ProvenanceStore,
) -> AsyncIterator[ExecuteEvent]:
    """Async generator: executes plan steps, yields ExecuteEvents."""
    # yield needed to make this an async generator
    for i, step in enumerate(plan.steps):
        started_at = datetime.datetime.utcnow().isoformat() + "Z"

        yield StepStartedEvent(
            step_index=i,
            primitive=step.primitive,
        )

        try:
            result: PrimitiveResult = await sandbox.run_primitive(
                step.primitive, step.params
            )
        except Exception as exc:
            logger.error(f"Step {i} ({step.primitive}) failed: {exc}")
            yield StepFailedEvent(
                step_index=i,
                primitive=step.primitive,
                error=str(exc),
            )
            # Record failure in provenance
            await store.record_step(
                session_id=plan.session_id,
                step_id=i,
                started_at=started_at,
                entry={
                    "primitive": step.primitive,
                    "params": step.params,
                    "input_hash": "",
                    "output_hash": "",
                    "duration_sec": 0.0,
                    "warnings": [],
                    "user_explanation": f"Step failed: {exc}",
                    "methods_paragraph": "",
                    "n_obs_before": None,
                    "n_obs_after": None,
                    "n_vars_before": None,
                    "n_vars_after": None,
                    "figures": [],
                },
            )
            # Halt on failure
            return

        # Record success in provenance
        warnings_dicts: list[dict[str, Any]] = [
            w.model_dump() for w in result.warnings
        ]
        await store.record_step(
            session_id=plan.session_id,
            step_id=i,
            started_at=started_at,
            entry={
                "primitive": result.primitive,
                "params": result.params,
                "input_hash": result.input_hash,
                "output_hash": result.output_hash,
                "duration_sec": result.duration_sec,
                "warnings": warnings_dicts,
                "user_explanation": result.user_explanation,
                "methods_paragraph": result.methods_paragraph,
                "n_obs_before": result.n_obs_before,
                "n_obs_after": result.n_obs_after,
                "n_vars_before": result.n_vars_before,
                "n_vars_after": result.n_vars_after,
                "figures": result.figures,
            },
        )

        yield StepCompletedEvent(
            step_index=i,
            result=result,
        )

        logger.info(
            f"Executor: step {i} ({step.primitive}) completed"
        )

    yield PlanCompletedEvent(n_steps=len(plan.steps))
    logger.info(f"Executor: plan completed ({len(plan.steps)} steps)")

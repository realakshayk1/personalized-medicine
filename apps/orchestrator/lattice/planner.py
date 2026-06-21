"""Planner: converts user messages to execution Plans via LLM.

The Planner is constrained to only emit PrimitiveInvocations for primitives
that are both in PRIMITIVE_REGISTRY and in the workflow's allowed_primitives.

For tests, inject a stub via `chat_fn` constructor argument.
"""

from __future__ import annotations

import json
import pathlib
from collections.abc import Callable
from typing import Any

import yaml
from loguru import logger
from pydantic import ValidationError

from lattice.models import Plan

# ---------------------------------------------------------------------------
# Workflow loader
# ---------------------------------------------------------------------------

_WORKFLOWS_DIR = pathlib.Path(__file__).parent.parent / "workflows"
_PROMPTS_DIR = pathlib.Path(__file__).parent.parent / "prompts"


def load_workflow(name: str) -> dict[str, Any]:
    """Load a workflow YAML by name."""
    path = _WORKFLOWS_DIR / f"{name}.yaml"
    if not path.exists():
        raise ValueError(f"Workflow '{name}' not found at {path}")
    with open(path) as f:
        return yaml.safe_load(f)  # type: ignore[no-any-return]


def load_system_prompt() -> str:
    """Load the planner system prompt markdown."""
    path = _PROMPTS_DIR / "planner_system.md"
    return path.read_text(encoding="utf-8")


def load_few_shot_examples() -> list[dict[str, Any]]:
    """Load few-shot planner examples from planner_examples.jsonl.

    Each line is a JSON object with keys: `user_message`, `anndata_summary`,
    and `plan` (a Plan-shaped dict). Blank lines are skipped. Returns an empty
    list if the file is absent so the planner still functions without examples.
    """
    path = _PROMPTS_DIR / "planner_examples.jsonl"
    if not path.exists():
        return []
    examples: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        examples.append(json.loads(line))
    return examples


def render_few_shot_block(examples: list[dict[str, Any]]) -> str:
    """Render few-shot examples as a prompt-appendable text block.

    Produces a deterministic, compact rendering of each example as a
    request/AnnData-summary -> Plan-JSON pair. Returns "" when there are no
    examples so the system prompt is unchanged in that case.
    """
    if not examples:
        return ""
    parts: list[str] = [
        "\n\n---\n\n## Few-shot examples\n",
        "Each example shows a user request (plus AnnData summary) and the exact "
        "JSON Plan you should emit. Follow these patterns; only propose "
        "primitives in the active workflow's allowed_primitives.\n",
    ]
    for i, ex in enumerate(examples, start=1):
        user_message = ex.get("user_message", "")
        summary = ex.get("anndata_summary", {})
        plan = ex.get("plan", {})
        parts.append(f"\n### Example {i}\n")
        parts.append(f"User request: {user_message}\n")
        parts.append(f"AnnData summary: {json.dumps(summary)}\n")
        parts.append(f"Plan:\n{json.dumps(plan)}\n")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Default chat function using Claude Agent SDK
# ---------------------------------------------------------------------------


async def _default_chat_fn(
    system_prompt: str,
    user_message: str,
) -> str:
    """Call Claude via the Claude Agent SDK and return the assistant text."""
    from claude_agent_sdk import query
    from claude_agent_sdk.types import (
        AssistantMessage,
        ClaudeAgentOptions,
        ResultMessage,
    )

    full_prompt = f"{system_prompt}\n\n---\n\nUser request: {user_message}"

    options = ClaudeAgentOptions(
        system_prompt="",  # system prompt embedded in prompt for simplicity
        tools=[],  # no tools needed — pure text completion
    )

    collected: list[str] = []
    async for msg in query(prompt=full_prompt, options=options):
        if isinstance(msg, AssistantMessage):
            for block in msg.content if hasattr(msg, "content") else []:
                if hasattr(block, "text"):
                    collected.append(block.text)
        elif isinstance(msg, ResultMessage):
            # ResultMessage carries the final text
            if hasattr(msg, "result") and msg.result:
                collected.append(str(msg.result))

    return "".join(collected)


# ---------------------------------------------------------------------------
# Planner
# ---------------------------------------------------------------------------


class UnknownPrimitiveError(ValueError):
    """Raised when the LLM emits a primitive name not in the registry/workflow."""


class Planner:
    """Converts user messages + AnnData summary + workflow into a validated Plan.

    Parameters
    ----------
    chat_fn:
        Async callable(system_prompt: str, user_message: str) -> str.
        Defaults to Claude Agent SDK. Pass a stub for tests.
    workflow_name:
        Name of the workflow YAML to use. Defaults to "standard_scrnaseq".
    """

    def __init__(
        self,
        chat_fn: Callable[..., Any] | None = None,
        workflow_name: str = "standard_scrnaseq",
    ) -> None:
        self._chat_fn: Callable[..., Any] = chat_fn or _default_chat_fn
        self._workflow_name = workflow_name
        self._workflow = load_workflow(workflow_name)
        self._few_shot_examples = load_few_shot_examples()
        # System prompt = primitive/ordering spec + appended few-shot examples.
        # The chat_fn seam receives this assembled prompt unchanged, so tests
        # that inject a stub chat_fn keep working.
        self._system_prompt = load_system_prompt() + render_few_shot_block(
            self._few_shot_examples
        )
        self._allowed = set(self._workflow.get("allowed_primitives", []))

    def _validate_plan(self, plan: Plan) -> None:
        """Raise UnknownPrimitiveError if any step references an unregistered primitive."""
        from lattice_primitives import PRIMITIVE_REGISTRY

        for step in plan.steps:
            if step.primitive not in PRIMITIVE_REGISTRY:
                raise UnknownPrimitiveError(
                    f"Primitive '{step.primitive}' is not registered in PRIMITIVE_REGISTRY. "
                    f"Registered: {list(PRIMITIVE_REGISTRY.keys())}"
                )
            if step.primitive not in self._allowed:
                raise UnknownPrimitiveError(
                    f"Primitive '{step.primitive}' is not in allowed_primitives for "
                    f"workflow '{self._workflow_name}'. Allowed: {sorted(self._allowed)}"
                )

    async def plan(
        self,
        session_id: str,
        user_message: str,
        anndata_summary: dict[str, Any],
    ) -> Plan:
        """Generate and validate a Plan from a user message.

        Parameters
        ----------
        session_id:
            The session to attach the plan to.
        user_message:
            Natural-language request from the user.
        anndata_summary:
            Dict with n_obs, n_vars, obs_columns, filename.

        Returns
        -------
        Plan — validated against PRIMITIVE_REGISTRY and workflow.
        """
        # Ensure primitives are registered
        import lattice_primitives.all_primitives  # noqa: F401

        context_msg = (
            f"Session ID: {session_id}\n"
            f"Workflow: {self._workflow_name}\n"
            f"AnnData summary: {json.dumps(anndata_summary)}\n\n"
            f"User request: {user_message}\n\n"
            f"Output ONLY a JSON Plan object (no markdown fences, no extra text)."
        )

        logger.info(f"Planner: generating plan for session {session_id}")
        raw_response = await self._chat_fn(self._system_prompt, context_msg)

        # Parse JSON — strip markdown fences if present
        cleaned = raw_response.strip()
        if cleaned.startswith("```"):
            lines = cleaned.split("\n")
            # Remove first and last fence lines
            cleaned = "\n".join(
                lines[1:-1] if lines[-1].strip() == "```" else lines[1:]
            )
        cleaned = cleaned.strip()

        try:
            plan_dict = json.loads(cleaned)
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Planner returned invalid JSON: {exc}\nRaw response: {raw_response[:500]}"
            ) from exc

        # Ensure session_id and workflow are set
        plan_dict["session_id"] = session_id
        plan_dict["workflow"] = plan_dict.get("workflow", self._workflow_name)

        try:
            plan = Plan.model_validate(plan_dict)
        except ValidationError as exc:
            raise ValueError(
                f"Planner output failed Pydantic validation: {exc}"
            ) from exc

        self._validate_plan(plan)
        logger.info(f"Planner: plan validated — {len(plan.steps)} step(s)")
        return plan

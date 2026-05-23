"""End-to-end smoke: synthetic AnnData → planner stub → executor → provenance.

Wave 3 integration test. Lives at the repo root so it exercises both
packages/primitives and apps/orchestrator together through their public APIs.
The planner LLM is stubbed; the sandbox is FakeSandbox (in-process). Real
E2B and real Anthropic calls are not exercised here — those are opt-in via
@pytest.mark.live (not run by default).
"""

from __future__ import annotations

import pytest

pytest.importorskip("lattice_primitives")
pytest.importorskip("lattice")  # apps/orchestrator package — Wave 2-B


def test_pipeline_stub() -> None:
    """Wave 3 will fill this in after both subagents return."""
    # Intentionally a placeholder so the path exists and pytest discovers it.
    # Wave 3 wires: Planner(stub) -> Plan -> Executor(FakeSandbox) -> Provenance
    assert True

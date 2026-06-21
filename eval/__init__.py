"""Lattice evaluation harness.

Three suites, mirroring AGENTS.md "Eval discipline":

- ``trap``  (instant): synthetic AnnData fixtures with injected failure modes.
  Asserts that the correct primitive sanity check fires. 100% recall is a hard
  gate — every trap MUST be caught or the suite fails.
- ``quick`` (~5 min): loads PBMC 3k, runs the standard scRNA-seq chain via
  ``run_primitive``, and computes a cluster-agreement placeholder metric against
  ``pbmc3k_processed`` labels when available. Network download failures SKIP
  gracefully rather than crashing.
- ``full``  (~90 min): currently runs quick + trap and documents the published
  datasets (5 treatment-vs-control studies, HCA immune subset) that are not yet
  wired in.

Run with::

    uv run python -m eval.run --suite quick
    uv run python -m eval.run --suite trap
    uv run python -m eval.run --suite full
"""

from __future__ import annotations

__all__ = ["run"]

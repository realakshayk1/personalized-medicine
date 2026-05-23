# Architecture Decision Records

ADRs are numbered, immutable once landed. To change a decision, write a new ADR that supersedes the old one. Reference each ADR in code and PRs that depend on it.

---

## ADR-001 — No free-form LLM code generation in the analysis path

**Date:** 2026-05-22
**Status:** Accepted (binding)

### Decision

The LLM may select primitives from the registry and fill their declared parameters. The LLM may **not**:

- Write the body of a primitive
- Modify the registry schema or any primitive's input/output schema
- Generate any code that is `exec`'d, `eval`'d, or otherwise executed against user data
- Bypass schema validation, sanity checks, or the provenance log
- Emit Scanpy/Seurat/R source code into the sandbox for execution

Every bioinformatics operation flows through `packages/primitives`. Primitives are human-written, version-pinned, schema-validated, and citation-bearing. The sandbox executes only registered primitives invoked with validated parameters.

### Reasoning

The canonical failure mode for LLM-driven biology analysis is the *silent case/control swap*: an LLM-generated Scanpy script that runs without error but produces inverted scientific conclusions, surfacing only at peer review. The constrained primitive registry is the architectural defense. See PRD.md §1.2 and BUILD_SPEC.md §3 Bet 1.

Secondary reasons: hallucinated Scanpy function names; missing `random_state`; unreproducible per-user code; intractable eval surface.

### Consequences

- Adding a new analysis capability requires adding a primitive (PR + tests + citation).
- The LLM cannot "be creative" — and that is the feature.
- Workflows that fall outside the registry are explicitly told "not supported yet; export your state and run manually."

### Enforcement

- Code review: any PR that introduces an `exec`/`eval`/`subprocess.run("python -c …")` path against user data is rejected without an ADR superseding this one.
- Test: `apps/orchestrator/tests/test_no_codegen.py` greps the orchestrator package for forbidden patterns and fails the build if any are present.
- The provenance log records every primitive call by registered name; any "anonymous" execution is impossible by construction.

---

## ADR-002 — Streaming protocol for execute endpoint

**Date:** 2026-05-22
**Status:** Proposed (subagent B finalizes)

The `/sessions/{id}/execute` endpoint streams per-step progress. SSE vs WebSocket is chosen by subagent B based on which maps more naturally onto the Claude Agent SDK's streaming output. Decision recorded here once made.

---

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
**Status:** Accepted (subagent B — SSE chosen)

### Decision

The `/sessions/{id}/execute` endpoint uses **Server-Sent Events (SSE)** via `sse-starlette`.

### Reasoning

1. The execute endpoint is strictly server→client (one-way). The client sends the plan once; the server streams events back. SSE maps onto this naturally; WebSocket's bidirectional channel adds complexity without benefit.
2. SSE maps onto the Claude Agent SDK's `async for message in query(...)` iterator directly: each iteration yields a message that becomes an SSE event.
3. SSE is HTTP/1.1 native, works through standard load balancers and proxies, and requires no special client-side library — the browser's `EventSource` API suffices.
4. SSE reconnect semantics (`Last-Event-ID`) give free retry without manual state management.
5. `sse-starlette` integrates trivially with FastAPI and adds ~200 lines of dependency vs. the complexity of `websockets`.

### Dependency

`sse-starlette>=2.1` added to `apps/orchestrator/pyproject.toml`.

---

## ADR-003 — Default sandbox is in-process (FakeSandbox); E2B is opt-in

**Date:** 2026-05-22
**Status:** Accepted (subagent B)

### Decision

The default sandbox for v0.1 local-first development is **FakeSandbox**, which executes primitives in-process via `lattice_primitives.run_primitive`. The `E2BSandbox` class exists but is never instantiated by default; it requires `E2B_API_KEY` and explicit opt-in.

### Reasoning

1. **Safety argument**: FakeSandbox is architecturally safe for v0.1 because all execution flows through the registered, schema-validated primitive registry (ADR-001). There is no `exec`/`eval` path. The only code that runs is human-written, version-pinned, citation-bearing primitive functions.
2. **Local-first**: v0.1 targets local desktop use (Tauri app). No network round-trip to E2B = lower latency, no cold-start cost, no API key required.
3. **HIPAA/data-residency**: user data never leaves the machine by default. Cloud burst is an explicit user choice.
4. **Simplicity**: FakeSandbox removes E2B cold-start latency (~5–15s), billing, and dependency on external service availability from the critical development path.
5. **Test isolation**: FakeSandbox is fully deterministic and requires no mocking of external services.

### When to use E2BSandbox

- Cloud burst: user requests processing on a remote server.
- Untrusted code execution (future: user-contributed primitives that haven't been audited).
- R primitives (DESeq2, sctransform) that require R 4.4 — not in v0.1 scope.

### Consequences

- All v0.1 tests use FakeSandbox (no E2B API key required in CI).
- Live E2B tests are gated with `@pytest.mark.live` and skipped unless `E2B_API_KEY` is set.
- The architectural guarantee (primitives = the only execution path) means the sandbox boundary is a performance and isolation boundary, not a security one, for v0.1.

---

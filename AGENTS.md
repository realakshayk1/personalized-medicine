# AGENTS.md — Lattice

> A README for coding agents. Read this first, every session.
>
> Lattice is an AI agent that lets bench biologists analyze scRNA-seq data through conversation. The product orchestrates a curated, version-pinned library of bioinformatics primitives — it never generates analysis code from scratch.

---

## The one invariant that matters

**Never let the LLM write analysis code. Ever.**

Every bioinformatics operation runs through a primitive in `packages/primitives`. Primitives have schemas, version pins, and sanity checks. The LLM picks **which** primitive to use and **what parameters** to fill — nothing else.

If you find yourself writing a code path where the model emits Scanpy/Seurat code that gets exec'd, stop and add a primitive instead. This is the architectural moat. See `docs/DECISIONS.md` ADR-001 for the reasoning.

---

## Project map

```
/lattice
  README.md               public 1-pager
  PRD.md                  product vision, users, GTM
  BUILD_SPEC.md           architecture, 12-week plan
  AGENTS.md               this file — read every session
  ROADMAP.md              dated, evolves weekly
  /docs
    PRIMITIVES.md         spec for every primitive (authoritative)
    EVAL.md               benchmark methodology + scores
    DECISIONS.md          ADRs for every meaningful tech choice
  /apps
    /desktop              Tauri v2 + Next.js 15 + TypeScript
    /orchestrator         FastAPI + Claude Agent SDK
  /packages
    /primitives           Python — the wrapped Scanpy/Seurat library
    /sdk                  TypeScript types shared frontend↔backend
  /eval
    /benchmarks           Tabula Muris, PBMC 3k, HCA immune, 5 papers
    /traps                datasets with known QC failure modes
```

---

## Tech stack (versions are not suggestions)

- **Python 3.12**, `scanpy>=1.11,<1.12`, `anndata>=0.10`, `scrublet`, `harmonypy`, `scvi-tools`, `celltypist`
- **R 4.4** inside E2B sandbox only — `DESeq2`, `Seurat 5`, `sctransform`. Called via `rpy2`. Never installed on user machine.
- **Node 20+**, **TypeScript 5.4+**, Next.js 15 App Router, **Tauri v2**, shadcn/ui, Tailwind v4
- **FastAPI 0.115+**, `claude-agent-sdk`, Pydantic v2, SQLAlchemy 2
- **E2B** (`e2b-code-interpreter` Python SDK) for sandboxed execution
- **Stripe** for billing, **SQLite** local + **Postgres 16** for cloud

Lock files: `uv.lock` for Python, `pnpm-lock.yaml` for JS. Don't `pip install` or `npm install` — use `uv add` and `pnpm add`.

---

## How to add a primitive

This is the most common task. Get it right.

A primitive lives in `packages/primitives/lattice_primitives/<category>/<name>.py`. Categories: `qc`, `preprocess`, `dim_reduce`, `cluster`, `integrate`, `annotate`, `de`, `plot`.

Every primitive needs:

1. The `@primitive` decorator with input/output schema, parameter spec, version pin, citation.
2. A function body that wraps **one Scanpy/Seurat call** (or a tightly-coupled pair like normalize_total + log1p). No business logic, no branching beyond the wrapped call.
3. A `@sanity_check` function that runs after execution and returns `List[Warning]`.
4. An explanation template (`templates/<name>.j2`) with two outputs: a one-sentence user message and a passive-voice Methods paragraph with citations.
5. Tests in `tests/primitives/test_<name>.py` covering: happy path, schema violation rejection, sanity check on known-bad input.

The reference implementation is `packages/primitives/lattice_primitives/preprocess/normalize_total_log1p.py`. Copy its shape exactly.

**Hard rules for primitives:**
- One Scanpy/Seurat call per primitive. If you need two operations, you need two primitives.
- Never mutate `.obs` columns that weren't declared in `produces`. Schema validator will reject this.
- Always set `random_state` from `params` if the underlying call accepts it. Reproducibility is non-negotiable.
- Never use `inplace=False` then assign — Scanpy semantics expect in-place by default. Follow the upstream convention.
- Citations go in the decorator, not in comments. The export pipeline reads them programmatically.

---

## How to add a new workflow (planner → executor)

A workflow is a *named sequence of primitives the planner is allowed to propose*. Workflows live in `apps/orchestrator/workflows/<name>.yaml`.

When adding one:
1. Define the workflow YAML with allowed primitives and ordering constraints.
2. Add few-shot examples to `apps/orchestrator/prompts/planner_examples.jsonl`.
3. Add eval cases to `eval/benchmarks/<workflow>/` with input AnnData and expected outputs.
4. Bump the planner prompt version in `apps/orchestrator/prompts/VERSION`.

Don't let the planner improvise outside declared workflows in v0.x. We'll loosen this later when evals justify it.

---

## Commands

```bash
# Setup
uv sync                              # Python deps
pnpm install                         # JS deps
pnpm --filter desktop dev            # Run Tauri app
uv run uvicorn lattice.app:app --reload  # Run orchestrator

# Tests — always run before commit
uv run pytest packages/primitives    # Primitive tests
uv run pytest apps/orchestrator      # Orchestrator tests
pnpm --filter desktop typecheck      # Frontend types
pnpm --filter desktop test           # Frontend tests

# Evals — run when touching primitives, planner, or models
uv run python -m eval.run --suite quick    # ~5 min, blocks PRs
uv run python -m eval.run --suite full     # ~90 min, run before releases

# Lint / format
uv run ruff check . --fix
uv run ruff format .
pnpm --filter desktop lint
```

After any commit that touches `packages/primitives` or `apps/orchestrator/prompts/`, run `uv run python -m eval.run --suite quick` and post results in the PR description. PRs without eval results are not mergeable.

---

## Code style — Python

- Type hints required, including for primitive parameters. `mypy --strict` in CI.
- Pydantic v2 for all schemas. No dataclasses for anything crossing a boundary.
- Loguru for logging, not the stdlib `logging` module.
- Function and variable names match Scanpy conventions where applicable (`adata`, `.obs`, `.var`, `.X`, `.uns`, `.obsm`). Don't rename them.
- AnnData mutations are explicit in the docstring; we follow Scanpy's in-place convention but document every column added.
- No `import *`. No `from x import *`.
- Async only at the FastAPI boundary; primitives are sync.

## Code style — TypeScript

- `strict: true` in tsconfig. No `any` without a `// reason:` comment.
- Server Components by default in Next.js; mark Client Components explicitly.
- shadcn/ui components live in `apps/desktop/components/ui` — extend, don't fork.
- Tailwind utility classes, not custom CSS. Design tokens in `tailwind.config.ts`.
- API routes are thin — they call into orchestrator over HTTP.
- Use the generated SDK types in `packages/sdk` rather than redefining shapes.

---

## What you decide vs. what you ask

**Decide on your own:**
- Implementation details inside an existing primitive
- File organization within established directories
- Test cases to add (more is better)
- Bug fixes that don't change public APIs
- Refactoring within a single file
- Docstring improvements

**Ask first:**
- Adding a new primitive — confirm category, schema, citation
- Changing a primitive's parameter schema (breaks reproducibility for existing sessions)
- Adding or removing a workflow
- Modifying provenance log structure
- Changing the planner prompt template
- Anything that touches billing, auth, or PII
- Anything that requires new dependencies

**Never do without explicit user instruction:**
- Run `pip install` or `npm install` directly — use `uv add` / `pnpm add`
- Modify version pins in `pyproject.toml` or `package.json`
- Push to main — feature branches only
- Force-push
- Delete or rewrite the provenance log schema
- Add free-form code generation paths to the orchestrator

---

## Failure modes to watch for

These are real and we have specific defenses against them. If you see code drifting toward any of them, push back.

1. **The case/control swap.** An LLM-generated script silently reverses treatment and control labels because the file naming was ambiguous. Schema validation on `.obs` columns is the defense. Never let a primitive overwrite a categorical condition column.

2. **The unpinned dependency.** A primitive works today, fails next month after a Scanpy minor release. Every primitive declares an exact version range. CI fails if `scanpy>=1.11,<1.12` is loosened without an ADR.

3. **The plausible hallucination.** Model emits `sc.tl.fancy_clustering()` which doesn't exist. The primitive registry is the allowlist — anything else is rejected.

4. **The dead-confidence cell type.** Single-LLM annotation says "memory CD8 T cell" with no markers cited. Cell type annotation requires consensus across two model calls + cited markers + Cell Ontology ID, or it's flagged "needs review."

5. **The 5GB context bomb.** Trying to send a 5M-cell AnnData into an LLM prompt. The orchestrator should never serialize raw matrices into prompts; always summary stats.

6. **The hidden randomness.** Forgot `random_state` somewhere; second run produces different UMAP. Every stochastic primitive accepts and threads a seed. Tests verify reproducibility.

---

## Eval discipline

Three benchmark suites in `/eval`:

- **quick** (~5 min): Tabula Muris 1K-cell subset, PBMC 3k. Runs on every PR.
- **full** (~90 min): all 5 published treatment-vs-control datasets, HCA immune subset, all trap datasets. Runs on release branches.
- **trap** (instant): synthetic datasets designed to trigger sanity-check failures. Lattice must surface the issue. 100% recall required.

Reported metrics live in `docs/EVAL.md` and are updated on every release. Regressions block merge.

Target numbers (do not loosen without an ADR):
- Cell-type annotation F1 vs published labels: ≥ 0.75
- DE top-50 overlap with published: ≥ 60%
- Trap detection recall: 100%
- Reproducibility (same input, same code → bit-identical output): 100%
- Standalone notebook validity (export → run on clean env → match): 100%

---

## Git workflow

- Branch per feature: `feat/<short-description>` or `fix/<short-description>`
- Conventional Commits: `feat(primitives):`, `fix(orchestrator):`, `chore(eval):`
- Squash on merge. Keep main linear.
- Open a PR even for solo work. Self-review the diff before requesting agent assistance.
- PR body must include: summary, what changed, what tested, eval-suite results if applicable.
- No commits to main directly. Ever.

---

## Working with the founders

This project is built by two undergraduates in parallel with school:
- **CS founder** (UF): owns frontend, Tauri, orchestrator, infrastructure
- **LSM founder** (Penn): owns primitives, eval, biology correctness, user research

Agent etiquette:
- When a change touches biology semantics (defaults, normalization choices, DE methods, annotation logic), surface it explicitly and tag the LSM founder for review. Don't just pick a "reasonable" default — biological defaults have downstream paper consequences.
- When a change touches Tauri, React, or infra, the CS founder owns the decision.
- For ambiguous changes that span both, ask. Don't pick.

---

## Pitfalls specific to this codebase

- **AnnData is mutable by default.** Scanpy's in-place convention means primitives that look pure aren't. Always document what gets added to `.obs` / `.var` / `.uns` / `.obsm` / `.layers`.
- **`adata.X` may be sparse or dense.** Many bugs come from assuming one. Use `scipy.sparse.issparse(adata.X)` checks where needed.
- **`rpy2` is not thread-safe.** R primitives run sequentially in the sandbox. Don't try to parallelize them.
- **E2B sandboxes have cold-start latency** (~5-15 sec cold, ~25ms warm). Keep one warm per session; pre-warm during file upload.
- **`var_names_make_unique()` should be called on ingest.** Duplicate gene symbols are common and break things silently downstream.
- **The user's `random_state` should propagate through the whole pipeline.** Set once in session, thread everywhere.
- **Matplotlib in a server context.** Always use `Agg` backend in the orchestrator. Never `plt.show()`.
- **Citations come from `packages/primitives/lattice_primitives/citations.yaml`**, not hardcoded strings. The export pipeline reads them programmatically.

---

## Reading order for new agents

If you have never seen this codebase before, read in this order before touching code:

1. This file (you're here)
2. `BUILD_SPEC.md` — architecture and the 12-week plan
3. `docs/PRIMITIVES.md` — the primitive contract
4. `docs/DECISIONS.md` — every ADR, especially ADR-001 (no free-form code gen)
5. One existing primitive end-to-end: `packages/primitives/lattice_primitives/preprocess/normalize_total_log1p.py` + its sanity check + its template + its test
6. `apps/orchestrator/workflows/standard_scrnaseq.yaml` — the canonical workflow

Then start with the smallest possible task. Add tests first; let the test be the spec.

---

*Last updated alongside BUILD_SPEC v0.1. When this file drifts from reality, fix this file first.*
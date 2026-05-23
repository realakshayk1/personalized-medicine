# Lattice — 90-second demo script

Day-14 milestone target. Browser-served (Tauri shell deferred).

## Setup (one-time)

```bash
uv sync
pnpm install
```

## Run

Terminal 1 — orchestrator:
```bash
uv run uvicorn lattice.app:app --reload --app-dir apps/orchestrator
```

Terminal 2 — frontend:
```bash
pnpm --filter @lattice/desktop dev
# Open http://localhost:3000
```

## Demo flow

1. Drag in a small `.h5ad` (PBMC 3k subset, or `eval/benchmarks/pbmc3k_500.h5ad` if downloaded).
2. Watch the inspection panel populate: cell count, gene count, obs columns.
3. Type: *"Run standard QC and preprocessing."*
4. Plan panel renders 4 steps:
   - calculate_qc_metrics
   - filter_cells_min_counts
   - filter_genes_min_cells
   - normalize_total_log1p
5. Click **Run Plan**. Steps tick green in order; sanity-check warnings render inline.
6. Open provenance: every step shows input hash, output hash, params, duration.

## What this proves (Day-14 acceptance)

- Frozen contracts (Pydantic + TS) survive end-to-end without modification.
- Constrained primitive registry holds: planner cannot invent primitives.
- Schema validation enforces case/control preservation through normalization.
- Provenance log captures full reproducibility envelope.

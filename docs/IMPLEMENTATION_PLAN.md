# Lattice — Implementation Plan: Everything Left for v0.1

> Status as of 2026-06-21. Target = BUILD_SPEC §1.3 "first workflow we ship": drag in `.h5ad` → conversational QC → clustering → annotation → DE → publication figures → exportable standalone notebook, validated by an eval suite hitting target metrics.
>
> This plan is grounded in current API research (scanpy 1.11.x, pydeseq2 0.5.x, decoupler, celltypist 1.7.x, mLLMCelltype, nbformat). API specifics are footnoted inline. Every new primitive spec below still needs LSM-founder sign-off per AGENTS.md before implementation.

---

## Where we are

**Done & solid:** primitive registry + schema validation + sanity-check framework (incl. case/control-swap guard), 4 primitives (`calculate_qc_metrics`, `filter_cells_min_counts`, `filter_genes_min_cells`, `normalize_total_log1p`), FastAPI orchestrator (planner→executor→provenance/SQLite, SSE), Next.js 4-panel UI, frozen SDK contracts, `FakeSandbox` (in-process), no-codegen AST test.

**The gap:** pipeline stops at *normalize*. 6 of 8 primitive categories empty (`dim_reduce`, `cluster`, `integrate`, `annotate`, `de`, `plot`). Zero eval infra. No notebook export. E2B is a stub. Planner prompt has no few-shot examples.

---

## Milestones (sequenced by demo value)

Each milestone is independently demoable. M1 is the highest-leverage step — it turns the QC demo into "drag a PBMC file in, get a labeled UMAP."

| # | Milestone | Unlocks | Rough size |
|---|-----------|---------|-----------|
| **M1** | Clustering chain: HVG → PCA → neighbors → UMAP → Leiden → UMAP plot | Labeled UMAP end-to-end | 5 primitives + 1 plot + figure rendering |
| **M2** | Eval harness + quick suite | PR gate, regression safety, AGENTS.md compliance | runner + 2 datasets + traps |
| **M3** | Marker DE + plots: `rank_genes_groups`, dotplot, heatmap, volcano | "what defines each cluster" | 1 DE + 3 plot primitives |
| **M4** | Notebook/report export | The reproducibility differentiator | export route + UI + .ipynb/.py |
| **M5** | Cell-type annotation: CellTypist + mLLMCelltype consensus | Named clusters with citations | 2 annotate primitives + CL mapping |
| **M6** | Cross-condition DE: pseudobulk + pydeseq2 | Treatment-vs-control | 1 aggregate + 1 DE primitive |
| **M7** | Integration: Harmony | Multi-batch datasets | 1 integrate primitive |
| **M8** | Real sandbox: wire E2B | Isolated execution (BUILD_SPEC bet) | sandbox completion |
| **M9** | Planner hardening | Correct primitive selection | few-shot examples + prompt |

Deferred to v0.2+ (do NOT build now): Tauri shell, Stripe/billing, auth, cloud/Postgres, collaboration, spatial, CRISPR.

---

## M1 — Clustering chain (highest priority)

New workflow `standard_scrnaseq` extends to: `calculate_qc_metrics → filter_cells_min_counts → filter_genes_min_cells → normalize_total_log1p → highly_variable_genes → pca → neighbors → umap → leiden → plot_umap`.

Each primitive copies the shape of `preprocess/normalize_total_log1p.py` exactly (decorator → one scanpy call → sanity_check → `templates/<name>.j2` → tests).

### 1.1 `dim_reduce/highly_variable_genes` ⚠️ biology-semantics
- **Call:** `sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key=None)`
- **requires:** `layer_state="log_normalized"` (for `seurat` flavor). *Note:* `seurat_v3` needs raw counts + `scikit-misc` — if we want v3 as default, HVG must run pre-normalization on a counts layer. **Decision needed: flavor + n_top_genes (paper consequences).** HVG is deterministic (no random_state).
- **produces:** `.var` cols `highly_variable`, `means`, `dispersions`, `dispersions_norm` (seurat); `layer_state` unchanged.
- **sanity:** flavor/data-state mismatch (seurat on raw counts, or v3 on log data); warn if `sum(highly_variable)` far from `n_top_genes`; assert `skmisc` importable if v3.
- **citation:** Stuart & Satija (Seurat v3) or Theis HVG — add to `citations.yaml`.

### 1.2 `dim_reduce/pca`
- **Call:** `sc.pp.pca(adata, n_comps=50, random_state=<seed>)` (uses HVGs automatically if `.var['highly_variable']` present — be explicit via `mask_var`).
- **requires:** `layer_state="log_normalized"`, `.var['highly_variable']` present.
- **produces:** `adds_obsm=["X_pca"]`, `.varm['PCs']`, `adds_uns=["pca"]`.
- **sanity:** NaNs in `X_pca`; warn if PC1 variance ratio dominant (>0.5).
- **random_state:** yes (default 0) — thread session seed.

### 1.3 `dim_reduce/neighbors`
- **Call:** `sc.pp.neighbors(adata, n_neighbors=15, n_pcs=None, random_state=<seed>)`.
- **requires:** `X_pca` in obsm.
- **produces:** `adds_uns=["neighbors"]`, `.obsp['distances']`, `.obsp['connectivities']`. *(Note: `.obsp` not in current AnnDataState schema — see "Schema gap" below.)*
- **sanity:** assert `X_pca` exists; warn if `n_neighbors >= n_obs`.

### 1.4 `dim_reduce/umap`
- **Call:** `sc.tl.umap(adata, random_state=<seed>)` (requires neighbors first; defaults `min_dist=0.5`, `init_pos='spectral'` are stable in 1.11).
- **requires:** `.uns['neighbors']`.
- **produces:** `adds_obsm=["X_umap"]`, `adds_uns=["umap"]`.
- **sanity:** NaNs in `X_umap`. *Reproducibility note:* numba multithreading can break bit-identical output even with a seed — eval may need single-threaded run.

### 1.5 `cluster/leiden` ⚠️ biology-semantics + version-critical
- **Call (future-proof, silences FutureWarning):** `sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False, random_state=<seed>)`. **Requires `python-igraph` dep.** Without explicit `flavor`, scanpy 1.11 warns the default will switch from `leidenalg` to `igraph`.
- **requires:** `.uns['neighbors']`.
- **produces:** `.obs['leiden']` (categorical), `adds_uns=["leiden"]`. **Decision needed: default resolution; BUILD_SPEC wants a multi-resolution sweep (0.2/0.5/0.8/1.2) — model as repeated calls with `key_added` or a separate sweep primitive.**
- **sanity:** assert neighbors exist; warn on under-clustering (one giant cluster) or over-clustering; confirm categorical.
- **dep risk:** some `python-igraph` versions throw int32 ValueError on large graphs (scanpy #3028) — pin & test on largest eval dataset.

### 1.6 `plot/plot_umap` + figure rendering pipeline
- **Call:** `sc.pl.umap(adata, color="leiden", show=False, return_fig=True)` → `fig.savefig(buf, format="png", dpi=150)` → `plt.close(fig)`.
- **Orchestrator:** set `matplotlib.use("Agg")` at import; never `plt.show()`.
- **produces:** nothing in AnnData; emits a figure. `PrimitiveResult.figures` already supports paths/data-URIs.
- **Frontend:** wire `FigurePanel.tsx` (currently placeholder) to render base64/URL figures streamed via SSE.
- **sanity:** assert buffer non-empty, figure has ≥1 axes.

### Schema gap to resolve first (touches frozen contract — needs ADR)
`AnnDataState` has `adds_layers/adds_obsm/adds_uns` but **no `adds_obsp` or `adds_varm`**. `neighbors` writes `.obsp`, `pca` writes `.varm['PCs']`. Options: (a) extend `AnnDataState` with `adds_obsp`/`adds_varm` (frozen-contract change → ADR-004), or (b) track these under `adds_uns` loosely. **Recommend (a) + ADR.** This is a prerequisite for M1.

---

## M2 — Eval harness (do early; it's the PR gate)

AGENTS.md mandates `uv run python -m eval.run --suite quick` on every PR touching primitives/prompts. The module doesn't exist. Build:

- `eval/run.py` (or `eval/__main__.py`) with `--suite {quick,full,trap}`, mirroring AGENTS.md.
- **quick suite:** `sc.datasets.pbmc3k()` (raw, for QC/preprocess) + `sc.datasets.pbmc3k_processed()` (has cluster→celltype labels for F1) + a Tabula Muris 1K subset.
- **trap suite (100% recall required):** synthetic AnnData with injected failures — high-mito, doublet contamination, batch confound, non-integer "raw" counts. Assert the right sanity-check fires.
- **Metrics:** annotation F1 ≥0.75, DE top-50 overlap ≥60%, trap recall 100%, reproducibility bit-identical, notebook validity 100%. Write results to `docs/EVAL.md`.
- **Datasets:** `scanpy.datasets` gives PBMC built-ins free. Tabula Muris from figshare (has `cell_ontology_class` + CL IDs). HCA immune via `cellxgene-census` (curated `cell_type_ontology_term_id`). Large data goes in `.gitignore` + a fetch script, not the repo.

---

## M3 — Marker DE + plots

- `de/rank_genes_groups` ⚠️ **Call:** `sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", tie_correct=True, pts=True)`. **Decision needed:** scanpy's default is `t-test`, not wilcoxon — best practice is wilcoxon; set explicitly (paper consequences). Writes `.uns['rank_genes_groups']`. Sanity: groupby categorical ≥2 groups, no size-1 group, log-normalized data.
- `plot/dotplot` — `sc.pl.dotplot(..., return_fig=True).savefig(...)` (returns composite DotPlot object, not Axes).
- `plot/heatmap` — `sc.pl.rank_genes_groups_heatmap` or `sc.pl.heatmap`.
- `plot/volcano` — **no built-in scanpy volcano**; build from `.uns['rank_genes_groups']` (`logfoldchanges` vs `-log10(pvals_adj)`) with raw matplotlib on an Agg figure.

---

## M4 — Notebook / report export (the differentiator)

- **Backend route** `/export/{session_id}` reads provenance log → builds notebook with `nbformat`: `new_notebook()`, per step a `new_markdown_cell` (explanation + Methods) + `new_code_cell` (the exact primitive call with locked params), header cell with input hash check, footer with citation list (read from `citations.yaml`). Set `kernelspec`.
- **`.py` export:** `jupytext --to py:percent` or `nbconvert --to python`.
- **Validity eval:** `nbconvert --execute` on a clean env, diff outputs (the "100% standalone validity" target).
- **Frontend:** `CodePanel.tsx` (stub today) shows the generated code; add export button.
- **New deps:** `nbformat`, `jupytext` (or `nbconvert`).

---

## M5 — Cell-type annotation ⚠️ heavy biology-semantics, LSM-founder owns

- `annotate/celltypist` — **Call:** `celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)`. **Hard input requirement:** log1p-normalized to exactly 1e4 (`normalize_total(target_sum=1e4)` then `log1p`) — different from our median default, so annotation needs its own normalized layer. New dep `celltypist`; models downloaded at runtime. Produces `.obs['predicted_labels']`, `.obs['majority_voting']`.
- `annotate/mllm_consensus` — multi-LLM consensus (mLLMCelltype pattern). Input = per-cluster top markers from `rank_genes_groups`. Opus + Sonnet (v0.1), independent annotation → consensus proportion (CP) + Shannon entropy → up to N deliberation rounds → CP>0.8 high / 0.5–0.8 review / <0.5 needs-expert. **This is the one sanctioned LLM-driven step** (BUILD_SPEC §2.5). Satisfies AGENTS.md failure-mode #4 (no dead-confidence cell types) only with cited markers + Cell Ontology IDs.
- **CL ID gap:** mLLMCelltype outputs free-text names, NOT `CL:XXXXXXX`. Need a separate name→CL mapper (CellOntologyMapper, or CL ontology lookup) to meet the "Cell Ontology ID or flagged needs-review" rule.
- This milestone is mostly LLM-orchestration logic, not a thin scanpy wrapper — design carefully so it still routes through a registered primitive (no free-form codegen).

---

## M6 — Cross-condition DE (pseudobulk + pydeseq2) ⚠️ biology-semantics

- `preprocess/pseudobulk` (aggregation) — `decoupler get_pseudobulk(adata, sample_col, groups_col, mode="sum", min_cells=10)` on **raw counts**. *(decoupler 1.x `dc.get_pseudobulk` vs ≥2.0 `dc.pp.pseudobulk`/`sc.get.aggregate` — pin and verify.)*
- `de/pseudobulk_deseq2` — **pydeseq2** (pure Python, avoids non-thread-safe rpy2): `DeseqDataSet(counts, metadata, design="~condition")` → `.deseq2()` → `DeseqStats(dds, contrast=["condition","treated","control"]).summary()`. **This is where the case/control-swap horror story lives — the `condition` immutable-column guard must be enforced and contrast orientation explicit.** Reserve R/DESeq2-in-E2B for parity checks only.
- New deps: `pydeseq2`, `decoupler`.

---

## M7 — Integration (Harmony)

- `integrate/harmony` — `sc.external.pp.harmony_integrate(adata, key="batch")` (harmonypy). requires `X_pca`; produces `.obsm['X_pca_harmony']`; downstream neighbors uses `use_rep="X_pca_harmony"`. Only triggered for multi-batch data (planner decides). harmonypy already in stack per AGENTS.md.

---

## M8 — Real E2B sandbox

- Complete `E2BSandbox.run_primitive` (error handling for malformed JSON output), wire a config flag to switch `FakeSandbox`↔`E2BSandbox`, pre-warm during upload, keep one warm per session (~25ms warm / 5–15s cold). Same execution path local & cloud. `FakeSandbox` stays the default/test path (ADR-003).

---

## M9 — Planner hardening

- Populate `apps/orchestrator/prompts/planner_examples.jsonl` with few-shot plans (standard QC+clustering, multi-batch, treatment-vs-control).
- Enrich `planner_system.md` with primitive descriptions + ordering constraints.
- Bump `prompts/VERSION`.
- Add eval cases that assert correct primitive selection on ambiguous prompts.

---

## Dependencies to add (needs explicit approval per AGENTS.md — `uv add`, never pip)

| Dep | For | Notes |
|-----|-----|-------|
| `scikit-misc` | HVG seurat_v3 | only if v3 chosen |
| `python-igraph` | leiden igraph flavor | pin; watch #3028 |
| `leidenalg` | leiden (current default) | one of igraph/leidenalg |
| `celltypist` | annotation | downloads models at runtime |
| `pydeseq2` | cross-condition DE | `design=` API (0.5.x) |
| `decoupler` | pseudobulk aggregation | verify 1.x vs 2.0 API |
| `nbformat` | notebook export | |
| `jupytext` or `nbconvert` | .py export / validity | |
| `cellxgene-census` | HCA eval data | eval-only, optional group |
| `harmonypy` | integration | already in AGENTS.md stack |

---

## Decisions that need a human (biology-semantics → LSM founder; contract → ADR)

1. **HVG flavor + n_top_genes** (seurat vs seurat_v3; 2000?) — paper consequences.
2. **Leiden default resolution** + whether to ship the multi-resolution sweep in v0.1.
3. **rank_genes_groups method** — wilcoxon (recommended) vs scanpy's t-test default.
4. **Annotation normalization** — CellTypist needs target_sum=1e4, conflicting with our median default; needs a dedicated layer.
5. **CL ID mapping approach** for annotation (which mapper).
6. **ADR-004: extend `AnnDataState`** with `adds_obsp`/`adds_varm` (frozen-contract change, blocks M1).
7. **igraph vs leidenalg** backend pin.

---

## Suggested order of execution

1. **ADR-004** (schema for obsp/varm) — unblocks everything in M1.
2. **M1** clustering chain + figure rendering — the demo moment.
3. **M2** eval harness — so M1 is regression-protected and PRs are mergeable.
4. **M3** marker DE + plots.
5. **M4** export — lock in the reproducibility story.
6. **M5/M6** annotation + cross-condition DE (the biology-heavy, LSM-owned milestones).
7. **M7/M8/M9** integration, real sandbox, planner polish.

After each milestone touching primitives: run `eval.run --suite quick`, post results in the PR (AGENTS.md gate).

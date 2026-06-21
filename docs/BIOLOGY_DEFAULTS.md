# Biology Defaults — Provisional (pending LSM-founder review)

> **Status: PROVISIONAL.** These are the default parameter choices baked into the
> v0.1 analysis primitives. Each has downstream paper consequences (AGENTS.md:
> "biological defaults have downstream paper consequences"). They are chosen from
> current best-practice research and are explicitly flagged for the LSM founder to
> confirm or override before we treat them as final. Changing a default after
> sessions exist breaks reproducibility, so finalize before wide use.
>
> Research basis: scanpy 1.11.x docs, Theis lab `sc-best-practices.org`,
> mLLMCelltype (bioRxiv 2025.04), the 2023–24 pseudobulk-DE consensus.

---

## How to read this

Each entry: **what we picked**, **why**, **the main alternative**, and **the
question for the LSM founder**. "Provisional default" = what ships now if nobody
overrides it.

---

## 1. HVG — highly variable genes (`dim_reduce/highly_variable_genes`)

- **Provisional default:** `flavor="seurat"`, `n_top_genes=2000`, `batch_key=None`.
- **Why:** `seurat` flavor operates on the **log-normalized** matrix we already
  produce (fits our `normalize → HVG` order with no extra step), is scanpy's
  default, and needs no extra dependency. 2000 is the field-standard count.
- **Main alternative:** `flavor="seurat_v3"` (the Seurat v3 VST method) is often
  preferred in recent best-practices, BUT it requires **raw counts** (so HVG must
  run *before* normalization on a counts layer) and the `scikit-misc` dependency.
  That's a larger pipeline change.
- **Question for LSM:** Is `seurat` on log-data acceptable for v0.1, or do we want
  `seurat_v3` on raw counts (reorders the pipeline + adds `skmisc`)? And is 2000
  the right count, or tissue-dependent?

## 2. PCA (`dim_reduce/pca`)

- **Provisional default:** `n_comps=50`, computed on HVGs only (`mask_var="highly_variable"`), `random_state` threaded from session seed (default 0).
- **Why:** 50 PCs is the scanpy/Theis convention; restricting to HVGs is standard
  and what downstream neighbors/UMAP expect.
- **Question for LSM:** 50 PCs default OK? Any datasets where we should surface an
  elbow/variance-ratio prompt instead of a fixed count?

## 3. Neighbors (`dim_reduce/neighbors`)

- **Provisional default:** `n_neighbors=15`, `n_pcs=None` (use all of `X_pca`),
  `metric="euclidean"`, `random_state=0`.
- **Why:** scanpy defaults; 15 is the canonical neighborhood size.
- **Question for LSM:** Defaults fine, or expose `n_pcs` selection tied to the PCA
  variance ratio?

## 4. UMAP (`dim_reduce/umap`)

- **Provisional default:** `min_dist=0.5`, `spread=1.0`, `init_pos="spectral"`,
  `random_state=0`.
- **Why:** scanpy defaults, stable through 1.11; visualization-only (does not feed
  clustering), so low biological risk.
- **Reproducibility note:** numba multithreading can make UMAP not bit-identical
  across machines even with a fixed seed. Eval may need single-threaded execution
  to hit the "bit-identical" target. Not a biology decision — flagged for eval.

## 5. Leiden clustering (`cluster/leiden`) — **highest-impact default**

- **Provisional default:** `resolution=1.0`, `flavor="igraph"`, `n_iterations=2`,
  `directed=False`, `random_state=0`. Single resolution (no sweep yet).
- **Why `igraph`:** scanpy 1.11 emits a FutureWarning unless you pass
  `flavor="igraph"`; the leidenalg default is being retired. `igraph` is faster and
  the documented future default. Requires the `igraph` package (added).
- **Why resolution 1.0:** scanpy default. BUT cluster count is highly sensitive to
  resolution, and BUILD_SPEC §1.3 step 8 calls for a **multi-resolution sweep
  (0.2 / 0.5 / 0.8 / 1.2) where the user picks**. We ship a single resolution now;
  the sweep is a deliberate follow-up.
- **Questions for LSM:**
  1. Default resolution 1.0, or lower (0.5 is common for cleaner top-level
     clusters)?
  2. Ship the multi-resolution sweep in v0.1, or single-resolution + manual
     re-run?
  3. `igraph` vs `leidenalg` backend — any reason to prefer leidenalg for parity
     with existing lab analyses?

## 6. Marker DE (`de/rank_genes_groups`) — **method choice matters**

- **Provisional default:** `method="wilcoxon"`, `tie_correct=True`, `pts=True`,
  `corr_method="benjamini-hochberg"`.
- **Why:** scanpy's *own* default is `t-test`, but the best-practices recommendation
  for marker detection is **Wilcoxon rank-sum** (non-parametric, robust). We set it
  explicitly rather than inheriting the t-test default.
- **Question for LSM:** Confirm Wilcoxon as the marker-DE default. (Note: this is
  *within-dataset marker detection*, distinct from cross-condition DE below.)

## 7. Cross-condition DE (`de/pseudobulk_deseq2`) — pseudobulk is non-negotiable

- **Provisional default:** aggregate to **pseudobulk** (sum counts per
  sample × cell type, `min_cells=10`) then **pydeseq2** (`design="~condition"`).
- **Why:** The 2023–24 single-cell community consensus is that naive per-cell
  Wilcoxon overstates significance (pseudo-replication); pseudobulk + a bulk method
  (DESeq2) is correct. pydeseq2 keeps it pure-Python and off the non-thread-safe
  rpy2 path; reserve R/DESeq2-in-E2B for parity checks.
- **Case/control guard:** the `condition` column is immutable-protected and the
  DESeq2 contrast orientation (`["condition","treated","control"]`) must be
  explicit — this is the exact failure mode AGENTS.md #1 defends against.
- **Questions for LSM:** Confirm pydeseq2 over rpy2/DESeq2 for v0.1. Default
  `min_cells` / minimum samples-per-group thresholds?

## 8. Cell-type annotation (`annotate/*`) — LSM owns end-to-end

- **CellTypist provisional:** `model="Immune_All_Low.pkl"`, `majority_voting=True`.
  Hard input requirement: log1p-normalized to **target_sum=1e4** (differs from our
  median normalization default → annotation needs its own normalized layer).
- **mLLMCelltype consensus provisional:** Opus + Sonnet, independent annotation →
  consensus proportion + Shannon entropy → ≤3 deliberation rounds; CP>0.8 high /
  0.5–0.8 review / <0.5 needs-expert. Cited markers required.
- **Open gap:** mLLMCelltype outputs free-text names, **not** Cell Ontology IDs
  (`CL:XXXXXXX`). AGENTS.md failure-mode #4 requires a CL ID or "needs review", so
  we need a name→CL mapping step (e.g. CellOntologyMapper / CL lookup).
- **Questions for LSM:** Default CellTypist model (immune vs tissue-specific)?
  Confidence thresholds for the consensus? Which CL-mapping approach?

## 9. Normalization (already shipped — for reference)

- **Current default:** `sc.pp.normalize_total(target_sum=None)` (median) + `log1p`.
  Confirmed still the Theis-lab shifted-log recommendation in 2025. No change
  proposed. Note the conflict with CellTypist's required 1e4 target (item 8).

---

## Summary table

| Step | Provisional default | Risk | LSM sign-off needed |
|------|---------------------|------|---------------------|
| HVG | seurat, 2000 genes | flavor choice (vs seurat_v3) | ✅ |
| PCA | 50 comps, HVG-only | low | optional |
| Neighbors | 15 neighbors, all PCs | low | optional |
| UMAP | min_dist 0.5, spectral | viz-only, low | no (eval: repro) |
| **Leiden** | **res 1.0, igraph, single** | **high (cluster count)** | ✅ + sweep decision |
| **Marker DE** | **wilcoxon** | **method choice** | ✅ |
| Cross-cond DE | pseudobulk + pydeseq2 | correctness-critical | ✅ |
| Annotation | CellTypist immune + consensus | LSM owns | ✅ |
| Normalization | median + log1p (shipped) | settled | no |

When these are confirmed, update this doc's status to FINAL and note any overrides
in the corresponding primitive's `ParamSpec` defaults + an ADR if a default changes
after sessions exist.

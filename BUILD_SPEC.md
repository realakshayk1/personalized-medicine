# Lattice — Build Specification

A concrete engineering plan for an AI agent that lets bench biologists run single-cell RNA-seq analysis through conversation.

This document covers **what to build and how to build it**, not how to pitch it. It answers the open questions from v0.3 with researched, opinionated decisions and gives you something you can start writing code against.

---

## 0. Decisions made up front (and the reasoning)

A list of every meaningful design decision in this spec, made before writing any code. If you disagree, change them here first, then change the spec.

| Decision | Choice | Why |
|---|---|---|
| Primary modality | scRNA-seq | Highest pain density; LSM co-founder is the user persona; standardized inputs (AnnData/h5ad); Theis lab `sc-best-practices` provides canonical reference workflow |
| File format | AnnData (`.h5ad`) as canonical, with auto-conversion from Seurat `.rds` and 10x Genomics raw output | AnnData has won; Scanpy is built around it; converters are mature |
| LLM provider | Claude Opus 4.7 for reasoning + Sonnet 4.6 for tool routing | Single-provider simplicity; Claude leads cell-typing benchmarks (bbaf677, Dec 2025) alongside GPT-5 and Kimi-k2; Agent SDK is built for this pattern |
| Agent framework | Claude Agent SDK (Python) | Anthropic's official package; handles tool loop, compaction, hooks; model-locked but that's fine because we're locked to Claude anyway |
| Code execution | E2B sandbox per session | 15M+ sandbox-sessions/month in production; "claude" template pre-built; ~25ms warm resume; isolates user data; lets us run arbitrary Scanpy/Seurat without touching our infra |
| Code generation strategy | Constrained primitive library, NOT free-form code gen | This is the trust differentiator. See §3. |
| Bioinformatics language | Python (Scanpy 1.11+) primary; R/Seurat via rpy2 for specific tools (DESeq2, sctransform) | Python ecosystem is broader; LSM co-founder likely knows Python better than R; AnnData native |
| Cell type annotation | LLM consensus (mLLMCelltype pattern, biorxiv April 2025) | 77.3% vs 61.3% for best traditional method on 50 datasets — this is the validated way |
| State storage | SQLite + AnnData files in user's local directory (Tauri app) OR S3 + Postgres (web app) | Start local-first with Tauri; ship web version once design partners want collaboration |
| Provenance | PROV-DM extension + Git-LFS for data hashes | Standard provenance model; biologists understand Git; reproducibility is the trust layer |
| Frontend | Next.js 15 + TypeScript + shadcn/ui | LSM co-founder can read it; UF CS partner can ship it fast |
| Desktop wrapper | Tauri v2 | Smaller than Electron; Rust core means file I/O doesn't choke on large h5ad files |
| Payment | Stripe + per-seat $29/mo personal, $99/mo lab seat (5+) | Below "ask PI for purchase" threshold; above "free" perception trap |
| Eval framework | Custom benchmark on Tabula Muris, Human Cell Atlas immune subsets, and 5 published treatment-vs-control datasets | Public ground truth; reproducible; published numbers we can beat |

---

## 1. The product, precisely

### 1.1 What a user does

The user is a bench biologist (grad student, postdoc, or junior staff scientist) who has just generated or downloaded an scRNA-seq dataset. They open Lattice. The flow:

1. Drag in an `.h5ad` file (or 10x Genomics output folder, or `.rds`).
2. Lattice auto-runs initial inspection: number of cells, number of genes, metadata columns, sample distribution. Surfaces this in a panel.
3. User types: *"Run standard QC and clustering. This is human PBMC, treated vs control."*
4. Lattice generates an **execution plan** (~6–12 steps), shows it as a checklist, and asks for confirmation before running.
5. On approval, Lattice executes step-by-step. Each step shows: the primitive invoked, its parameters, the code run, the output (figure/table), and a sanity-check result.
6. At the end, user gets a notebook view: full analysis, all figures, downloadable as a standalone `.ipynb` or `.py` script that runs without Lattice.
7. User can ask follow-up questions: *"Why is cluster 4 mostly from the treated sample?"* or *"Re-cluster with leiden resolution 0.5 and tell me what changed."*
8. The session persists. Reopen tomorrow, everything is exactly where you left it.

### 1.2 What it should NEVER do

These are firm rules. Violations are bugs.

- **Never write analysis code from scratch.** Lattice orchestrates primitives. The primitives are version-pinned and audited.
- **Never silently swap case/control labels.** Schema validation prevents the canonical horror story.
- **Never hide what it did.** The code is always visible, downloadable, and reproducible standalone.
- **Never invent biology.** Cell type annotations cite their source markers; pathway claims cite databases; gene functions cite NCBI Gene / UniProt.
- **Never run analyses that take >10 min without explicit user approval** (with cost/time estimate beforehand).
- **Never upload user data without consent.** Local-first by default. Cloud is opt-in.

### 1.3 The first workflow we ship

End-to-end scRNA-seq from `.h5ad` to publication figures. The specific steps Lattice can execute in v0.1:

1. **Ingest** — load h5ad, validate, summarize
2. **QC** — flag mitochondrial-heavy cells, doublets (Scrublet), low-count cells; **always plot before filtering, then plot after**
3. **Filter** — apply user-confirmed thresholds (defaults derived from data, not hardcoded)
4. **Normalize** — shifted log transform (Theis lab recommendation), with option for scran or sctransform
5. **HVG selection** — `seurat_v3` flavor, top 2,000 by default
6. **Dimensionality reduction** — PCA → neighbors → UMAP
7. **Integration** (if multi-batch) — Harmony default; offer scVI for complex cases
8. **Clustering** — Leiden with multi-resolution sweep (0.2, 0.5, 0.8, 1.2); user picks
9. **Cell type annotation** — CellTypist + LLM consensus (mLLMCelltype-style) with citations
10. **Differential expression** — `rank_genes_groups` (Wilcoxon by default; t-test or logreg on request); cross-condition DE via `pseudobulk + DESeq2` when `.obs` has sample IDs (the correct approach per 2024 single-cell community consensus)
11. **Figures** — publication-defaults: UMAP, marker dotplot, volcano plot, heatmap of top DE genes
12. **Report** — auto-generated Markdown with figures, tables, code, and a "Methods" section the user can paste into a paper

Every step has schema in, schema out, sanity-check outputs, and a one-liner explanation surfaced to the user.

---

## 2. Architecture

### 2.1 Components

```
┌────────────────────────────────────────────────────────────────┐
│  TAURI DESKTOP APP (Rust shell + Next.js frontend)             │
│  • File picker for .h5ad / 10x / .rds                          │
│  • Chat panel + plan panel + figure panel + code panel         │
│  • Local SQLite for sessions, settings, billing token          │
└────────────────────────────────────────────────────────────────┘
                  │
                  │ (gRPC or local HTTP)
                  ▼
┌────────────────────────────────────────────────────────────────┐
│  ORCHESTRATOR SERVICE (Python, FastAPI)                        │
│  Built on Claude Agent SDK                                      │
│  • Planner: parses user intent → execution plan                │
│  • Executor: walks the plan, invokes primitives                 │
│  • Verifier: runs sanity checks after each step                 │
│  • Reporter: streams progress to frontend                       │
└────────────────────────────────────────────────────────────────┘
                  │
                  ├─────────────┐
                  ▼             ▼
┌──────────────────────┐  ┌─────────────────────────────────────┐
│ PRIMITIVE REGISTRY    │  │ EXTERNAL KNOWLEDGE                  │
│ (Python package)      │  │ • PanglaoDB markers                  │
│ Version-pinned tools  │  │ • CellTypist model registry          │
│ • Scanpy 1.11.x       │  │ • Cell Ontology (CL)                 │
│ • Scrublet 0.2.x      │  │ • MSigDB pathways                    │
│ • Harmony, scVI       │  │ • GSEApy enrichment                  │
│ • DESeq2 (via rpy2)   │  │ • LiteratureRAG: PubMed via Entrez   │
│ • CellTypist          │  └─────────────────────────────────────┘
│ • Scanpy plotting     │
└──────────────────────┘
                  │
                  ▼
┌────────────────────────────────────────────────────────────────┐
│  E2B SANDBOX (one per session)                                 │
│  • Python 3.12, scanpy 1.11, scrublet, harmony, scVI, rpy2     │
│  • R 4.4 with DESeq2, Seurat 5                                 │
│  • User's AnnData mounted at /workspace/data.h5ad              │
│  • All execution happens here, isolated                         │
└────────────────────────────────────────────────────────────────┘
                  │
                  ▼
┌────────────────────────────────────────────────────────────────┐
│  PROVENANCE STORE (SQLite local; Postgres if cloud)            │
│  • Every primitive invocation logged                           │
│  • Input hash, output hash, parameters, version, timing         │
│  • Auto-generates standalone notebook on export                 │
└────────────────────────────────────────────────────────────────┘
```

### 2.2 Why this stack specifically

**Tauri over Electron**: 600MB Electron app vs 20MB Tauri app. Scientists with 50GB h5ad files care about disk and RAM headroom.

**Claude Agent SDK over LangChain/LangGraph**: less abstraction tax, built-in tool loop, hooks for guardrails, model-locked to Claude (a feature). The 2026 production guides repeatedly note that the Agent SDK is faster to ship and harder to misuse than the alternatives.

**E2B over direct subprocess**: industry standard for AI-generated code execution. The "claude" template at E2B is already configured. Microsecond startup (well, ~25ms warm resume). When your agent occasionally writes broken code, the user's machine doesn't crash — the sandbox does. Critical for trust: a single `rm -rf` incident (this has happened with Claude Code users) is reputation-killing.

**Scanpy primary over Seurat**: the Theis lab `sc-best-practices` book is essentially the canonical reference for the modern field, and it's Scanpy-first. Seurat support added via rpy2 for specific tools where R is dominant (DESeq2, sctransform, edgeR).

**SQLite + local files first**: undergrad founders shouldn't be running production cloud infrastructure. Local-first means privacy is the default, no HIPAA, no GDPR issues, and the architecture is dead simple. Cloud collaboration is a v2 feature.

### 2.3 The primitive registry — the hardest design decision

This is the technical IP. Read carefully.

A **primitive** is a wrapper around a single Scanpy or Seurat function with:

1. **Strict input schema**: required `.obs` columns, required dtype, required state of the AnnData (e.g., "must have run `normalize_total` before this")
2. **Strict output schema**: what gets added to `.obs`, `.var`, `.uns`, `.obsm`
3. **Parameter schema**: types, ranges, defaults, descriptions
4. **Sanity check function**: runs after execution, returns warnings if outputs look wrong
5. **Explanation template**: a Jinja string that surfaces to the user what was done
6. **Version pin**: explicit Scanpy version, model version, parameter version

Example pseudocode:

```python
@primitive(
    name="normalize_total_log1p",
    requires=AnnDataState(
        obs_columns=["total_counts"],  # must have run calculate_qc_metrics
        layer_state="raw_counts"
    ),
    produces=AnnDataState(
        layer_state="log_normalized",
        adds_layers=["counts_normalized"]
    ),
    params={
        "target_sum": ParamSpec(
            type=float,
            default=None,  # None = use median; this is the Theis lab default
            description="Per-cell scaling target. None uses median of total_counts (recommended). 1e4 = CP10k, 1e6 = CPM."
        ),
        "exclude_highly_expressed": ParamSpec(
            type=bool,
            default=False,
            description="Exclude top-1% genes from normalization. Use for tissue with dominant cell types."
        )
    },
    scanpy_version=">=1.11.0,<1.12",
    citation="Luecken & Theis, Mol Sys Biol 2019"
)
def normalize_total_log1p(adata: AnnData, target_sum=None, exclude_highly_expressed=False) -> AnnData:
    adata.layers["counts_normalized"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed)
    sc.pp.log1p(adata)
    return adata

@sanity_check("normalize_total_log1p")
def check_normalization(adata_before, adata_after, params) -> List[Warning]:
    warnings = []
    # Did the median total_counts roughly equal target_sum?
    median_after = np.median(adata_after.X.sum(axis=1))
    expected = params["target_sum"] or np.median(adata_before.obs["total_counts"])
    if abs(median_after - np.log1p(expected)) > 1.0:
        warnings.append(Warning(
            severity="medium",
            message=f"Median normalized total ({median_after:.2f}) deviates from expected ({np.log1p(expected):.2f})"
        ))
    # Is the data actually log-transformed? Max should be < ~10
    if adata_after.X.max() > 15:
        warnings.append(Warning(
            severity="high",
            message="Maximum value after log1p exceeds 15 — data may not be raw counts"
        ))
    return warnings
```

The LLM never writes the body of `normalize_total_log1p`. The LLM picks **which primitive to use** and **what parameters to fill**. The body is human-written, version-pinned, schema-validated.

This is what kills the "LLM silently reversed case/control" failure mode. The schema enforces that the `condition` column in `.obs` is preserved through every operation. Any primitive that touches `.obs` declares which columns it modifies.

### 2.4 The "explain what you did" layer

Every primitive emits, alongside its output:

- **One sentence for the user**: "Filtered out 1,247 cells with high mitochondrial content (>15%) and 89 likely doublets (Scrublet score >0.25). 14,302 cells remain."
- **Methods-section paragraph for paper**: written in passive scientific voice with parameter values and citations.
- **Sanity check summary**: green checkmarks for things that look fine, yellow warnings for things to consider, red flags for likely problems.

The LLM **synthesizes these from templates**, not from scratch. The template guarantees the parameter values and citations are correct.

### 2.5 Cell type annotation — the mLLMCelltype pattern

For step 9 (cell type annotation), Lattice uses a **multi-LLM consensus** following the mLLMCelltype framework (bioRxiv April 2025, 77.3% mean accuracy across 50 datasets):

1. After clustering, for each cluster, extract top-N differential markers via `rank_genes_groups`.
2. Build a structured prompt: "Cluster X has top markers: [gene1, gene2, ...]. Tissue context: [user-provided]. Possible cell types from this tissue: [from PanglaoDB filtered by tissue]."
3. Send to Claude Opus and Claude Sonnet (different temperatures) in parallel.
4. If both agree, return the consensus annotation with confidence "high."
5. If they disagree, run a deliberation round: each model sees the other's annotation and reasoning, can revise.
6. After 2 rounds, if still divergent, flag for user review with the disagreement explicit.
7. Always output the supporting markers, the cell ontology ID (CL:XXXXXXX), and the source paper for the canonical marker set.

This is one of the few places Lattice does **LLM-driven annotation** rather than pure tool orchestration — because LLMs measurably beat traditional methods here (~15% accuracy advantage). The trust comes from the **transparent reasoning chain and consensus**, not from blind LLM output.

### 2.6 The reproducibility layer

Every session, in the background, accumulates a structured execution log:

```json
{
  "session_id": "uuid",
  "started": "2026-05-22T14:00:00Z",
  "input_data": {
    "path": "user/data.h5ad",
    "sha256": "abc123...",
    "n_obs": 16934,
    "n_vars": 36601
  },
  "steps": [
    {
      "step_id": 1,
      "primitive": "calculate_qc_metrics",
      "version": "scanpy==1.11.1",
      "params": {"percent_top": [50, 100, 200, 500], "log1p": false},
      "input_hash": "abc123...",
      "output_hash": "def456...",
      "duration_sec": 2.4,
      "warnings": [],
      "user_explanation": "Computed QC metrics for each cell."
    },
    // ...
  ]
}
```

On "Export Notebook", Lattice generates a standalone `.ipynb`:

- Cell 1: data load with hash check
- Cells 2-N: each primitive call as a code cell, with markdown above explaining what it does
- Cell N+1: Methods-section paragraph
- Cell N+2: citation list

The user runs the notebook outside Lattice and gets bit-for-bit identical results.

**This is the differentiator from generic AI code generation**: a biologist can hand the notebook to a reviewer or PI and the reviewer can re-run it. The Substack horror story (LLM silently reversed case/control, paper went to peer review with wrong conclusions) literally cannot happen if every parameter, every version, every input hash is locked.

---

## 3. The two non-obvious technical bets

### Bet 1: Constrained primitives, not free-form code generation

The instinct is "let the LLM write Scanpy code." Don't.

**Why it loses**:
- LLMs hallucinate Scanpy function names that don't exist (this happens regularly with even Claude/GPT-4)
- LLMs forget to set `random_state`, breaking reproducibility
- LLMs invert axes, swap labels, and silently corrupt data — and the resulting code *runs without error*
- Every user gets slightly different code; no eval surface
- No way to verify correctness short of reading every generated script

**Why constrained primitives win**:
- Every analysis goes through a finite, audited set of functions
- Eval is tractable: 50 primitives × 10 test datasets × 5 parameter combinations = a real benchmark suite
- Bugs are fixed once for all users
- When Scanpy updates, you update primitives in one place
- The LLM does what it's good at (parameter selection, plan construction, explanation) and doesn't do what it's bad at (writing novel scientific code from scratch)

**The tradeoff**: you can't support arbitrary user workflows in v0.1. Someone wanting a weird custom analysis is told "Lattice doesn't support that yet — here's how to export your data and run it yourself."

### Bet 2: Local-first, sandboxed execution per session

The instinct is "build a cloud SaaS so you can scale and track usage."

**Why local-first wins for v0.1**:
- No HIPAA, no GDPR, no SOC 2 needed to ship — your users' data never leaves their machine
- PIs love this — they don't have to argue with their IRB or institution
- Tauri + E2B local mode means the same code runs on the user's laptop or in cloud
- Eliminates a class of "your servers are down" customer issues
- Two undergrads can't run reliable cloud infrastructure; don't try
- Local-first is the right answer for academic data anyway

**Why E2B for the execution sandbox (even locally)**:
- Even on a user's machine, you want the AI-generated/-selected code in an isolated process
- E2B local mode uses gVisor; the user's machine is protected from `rm -rf` accidents
- Same execution path for local and cloud means you don't fork the code

**When you go cloud**: when a customer asks. v2 feature, not v0.1.

---

## 4. The 12-week build plan, day-by-day-ish

### Week 1: Skeleton
- **Day 1-2**: monorepo (Turborepo); `apps/desktop` (Tauri + Next.js); `apps/orchestrator` (FastAPI); `packages/primitives` (Python); `packages/sdk` (TypeScript types shared)
- **Day 3-4**: Tauri app loads, file picker accepts .h5ad, displays AnnData summary
- **Day 5**: Orchestrator runs Claude Agent SDK with one trivial tool ("describe_anndata")
- **Day 6-7**: end-to-end: drag h5ad → orchestrator reads it → Claude sees summary → returns text → frontend displays

### Week 2: First three primitives + E2B
- **Day 8-9**: E2B sandbox integration; orchestrator launches a sandbox per session; AnnData uploaded into sandbox
- **Day 10-11**: implement primitive registry (decorator, schema, sanity check); first three primitives: `calculate_qc_metrics`, `filter_cells_basic`, `normalize_total_log1p`
- **Day 12**: planner agent: takes user prompt, returns a plan as an ordered list of primitive invocations
- **Day 13-14**: executor: walks the plan, calls each primitive in the sandbox, streams progress

### Week 3: Clustering and visualization
- **Day 15-16**: 7 more primitives: `highly_variable_genes`, `pca`, `neighbors`, `umap`, `leiden_clustering`, `rank_genes_groups`, `plot_umap`
- **Day 17-18**: figure rendering: matplotlib output → PNG/SVG → frontend display
- **Day 19**: verifier agent: runs sanity checks after each step, surfaces warnings
- **Day 20-21**: end-to-end test on Tabula Muris subset — drag in, ask for "standard clustering", get UMAP

### Week 4: Cell type annotation
- **Day 22-23**: PanglaoDB ingestion; CellTypist integration
- **Day 24-25**: mLLMCelltype-style consensus annotation (Opus + Sonnet, 2 rounds max)
- **Day 26-27**: citation rendering — every cell type label links to its source markers and Cell Ontology ID
- **Day 28**: end-to-end: PBMC dataset → 8 named clusters with citations

### Week 5: Differential expression and pseudobulk
- **Day 29-30**: `rank_genes_groups` primitive with multiple methods; volcano plot primitive
- **Day 31-33**: pseudobulk + DESeq2 (via rpy2) — this is the modern correct way to do cross-condition DE
- **Day 34-35**: cross-condition DE workflow: cluster, pseudobulk per sample, DESeq2

### Week 6: Polish and design partners
- **Day 36-37**: notebook export — generate runnable `.ipynb` with full provenance
- **Day 38-39**: report export — auto-generated Markdown with figures, tables, methods paragraph
- **Day 40-42**: onboard first 3 design partners from Penn LSM network; sit in the room while they use it

### Week 7-8: Listen and fix
- Whatever the design partners hit, fix
- Add: bulk RNA-seq support (DESeq2 from count matrix), Harmony integration for batch effects
- Add: doublet detection (Scrublet)

### Week 9-10: Trust layer hardening
- Provenance UI: show the full execution log; let users diff between two sessions
- Sanity-check coverage: every primitive has 3+ checks
- The "Methods" section export becomes paper-grade

### Week 11: Pricing and Stripe
- $29/mo personal, $99/seat/mo lab tier (5+ seats)
- Free tier: 5 sessions/month, no export

### Week 12: Demo-ready
- 25+ active design partner users at Penn + UF
- 5+ paying ($25-99/mo)
- 1 testimonial captured
- A 90-second demo video where someone analyzes a real dataset start to finish

---

## 5. Eval methodology

This is what makes the difference between a toy and something a real lab will pay for.

### 5.1 The benchmark dataset

Build a private repository with:

1. **Tabula Muris (mouse cell atlas)** — published cell types, 100K cells, multiple tissues. The "easy" baseline.
2. **PBMC 3k (10x Genomics)** — published annotation; the most-replicated dataset in the field; minimal floor.
3. **Human Cell Atlas immune subset** — fine-grained T-cell subtypes; tests cell-type annotation precision.
4. **5 published treatment-vs-control datasets** from recent papers (immune checkpoint, CAR-T, drug perturbation). Source the AnnData + the supplementary tables. Compare Lattice's DE results to published DE.
5. **2 "trap" datasets** — datasets with known QC failure modes (high mito, doublet contamination, batch effects). Lattice should *flag these* not silently process them.

### 5.2 What to measure

For each dataset, after Lattice runs end-to-end without user intervention:

- **Cell type annotation accuracy**: F1 vs published labels at cluster level. Target ≥0.75 (mLLMCelltype gets 0.77 across 50 datasets; we should match or beat with primitive consistency).
- **DE overlap with published**: top-50 DE genes per comparison, what fraction overlap with the published paper. Target ≥60% in top-50.
- **QC flagging recall**: on the trap datasets, did Lattice surface the known issue? Target 100% — these issues are not subtle.
- **Reproducibility**: run the same dataset twice. Bit-for-bit identical outputs? With `random_state` pinned, yes.
- **Standalone notebook validity**: export the notebook, run it on a fresh machine, do outputs match? Target 100%.

### 5.3 What to publish

By week 12, publish a preprint comparing Lattice on these benchmarks vs:
- Manual Scanpy analysis (by an experienced bioinformatician — your LSM co-founder)
- Generic ChatGPT/Claude with no scaffolding ("here's an h5ad file, write me a clustering script")
- Galaxy single-cell workflow

This becomes the credibility currency. "We are accurate AND faster AND reproducible" beats "we are cool AI" any day.

---

## 6. The most likely failure modes

| Failure | When it happens | Mitigation |
|---|---|---|
| LLM picks the wrong primitive ("use this when you should have used that") | When user prompts are ambiguous, or when the planner doesn't have enough context | Verifier agent reads the output and flags weirdness; user can revert any step; planner asks clarifying questions for ambiguous inputs |
| Sanity check is too noisy → users ignore warnings | Always-present yellow warnings train users to dismiss | Three-tier severity (info/warn/error); only errors block; warnings shown once per session |
| User wants to do something not in the primitive library | Every week | Two options surfaced: (a) "Here's how to do this manually, exported as a notebook starting from your current state" (b) waitlist for a custom primitive — feedback channel directly to engineering |
| The h5ad file is malformed or 50GB | All the time, surprisingly | Validate on ingest; refuse to load if invalid; for huge files, prompt for cloud burst mode |
| Cell type annotation is confidently wrong | When the cluster has weird markers (ambient RNA, dying cells) | Consensus mechanism flags low-confidence clusters; "unannotated / needs expert review" is a valid final state |
| Sandbox cold start is 20 seconds and users hate it | Every new session | Keep one warm sandbox per user; pre-warm during file upload |
| The user's R installation breaks rpy2 calls | All the time, R is fragile | R only runs inside E2B sandbox where we control the environment; user never installs R |
| LLM token cost per analysis is $5 and margin is bad | At scale | Use Sonnet for routing/parameter selection (~$3/M tokens); only use Opus for cell-type consensus and synthesis; budget cap per session with user warning |
| Two design partners want opposite features | Week 6+ | Pick one; document the trade-off; the other becomes v0.2 |

---

## 7. What you should personally validate this week

Six concrete things to do this week, in order:

1. **Run Tabula Muris through Claude Opus manually with the Claude Agent SDK, no product wrapper, no primitives** — just give Claude the Scanpy docs and an h5ad file and see what it does. Where does it fail? This is your real Day 0 calibration. Three hours.

2. **Read the Theis lab `sc-best-practices` book** (sc-best-practices.org) cover-to-cover. This is the reference manual for the field. Your primitive library should match it. Two evenings.

3. **Install E2B and run their "claude" template** end-to-end. Make sure you can launch a sandbox, upload a file, run Scanpy in it, and get results back. One afternoon.

4. **Pick 5 published scRNA-seq papers from the past 12 months** in immunology/cancer/neuro. Download their data from GEO. These are your benchmark gold standards. One afternoon.

5. **Interview 5 biology grad students at Penn** — your LSM co-founder's classmates — about their scRNA-seq workflow. Specifically: how long does it take them, where do they get stuck, who do they ask for help, what do they wish existed. Three hours total.

6. **Sketch the primitive registry decorator** in Python. Implement 3 primitives end-to-end against a test h5ad. Make sure the schema validation actually catches a case/control label swap. One day.

If any of these six surfaces something that breaks the plan, change the plan before writing the rest of the code.

---

## 8. Open questions I'm flagging but not resolving

These are real questions that I don't have enough information to decide on. Bring them to a v0.4 conversation:

- **R/Seurat support depth**: how much do we lean on rpy2 for R tools vs. building Python-only with eventual native ports? Probably depends on the design-partner mix. Some labs are 100% R.
- **Multi-modal**: when does spatial transcriptomics (Visium, Xenium) become important enough to add? Adoption is real but the field is still settling on standard formats.
- **Long-context for very large datasets**: a 5M-cell dataset has too many cells to discuss "by cluster" with the LLM. Architectural question: how do you let the LLM reason about millions of cells without dumping them all in context? Probably summary statistics + targeted sampling, but it's a real design problem.
- **Collaboration features**: when does "lab" mean "shared workspace with permissions" vs. "we all bought seats"? v2 question.
- **CRISPR screen analysis**: MAGeCK is the standard but pre-LLM. Adding it as a workflow could open up a different buyer (functional genomics labs). Not urgent.
- **Bulk RNA-seq UI**: the same orchestrator + primitive pattern works, but the workflow is meaningfully different (count matrix in, DESeq2 out). Worth thinking about as a separate sub-product or as a workflow within Lattice.

---

## 9. The pre-flight checklist before the first commit

Before writing line one of code, you should be able to answer yes to all of these:

- [ ] Tabula Muris loads in Scanpy locally without error on both your machines
- [ ] You've manually clustered a PBMC dataset in Scanpy and you both understand each step
- [ ] You have an Anthropic API key with sufficient credits ($200+) for development
- [ ] You have an E2B account and have launched a sandbox at least once
- [ ] You've read the Theis lab QC and normalization chapters at minimum
- [ ] You've talked to at least 5 biology grad students about their scRNA-seq pain
- [ ] You've decided who owns frontend (Next.js + Tauri) and who owns orchestrator (Python + Agent SDK)
- [ ] You've agreed on a Git workflow and weekly cadence given school schedules

Day 1 of coding starts after this list is complete.

---

## Appendix: Decision points where the research disagrees with my v0.3 PRD

The research surfaced things that change earlier recommendations:

1. **Cell type annotation should use LLM consensus, not just CellTypist.** mLLMCelltype framework, bioRxiv April 2025, validated 77.3% mean accuracy across 50 datasets, 8M cells. Beats traditional methods by ~15%. This is one of the few places LLMs are demonstrably better than dedicated tools — embrace it.

2. **Multi-LLM (not just single Claude) wins for cell-type annotation.** The benchmark (Briefings in Bioinformatics, Dec 2025) shows that Kimi-k2, GPT-5, Claude-4.1, and Grok-4 all perform similarly well — and consensus across them is better than any single one. For v0.1 use Opus + Sonnet (both Claude); for v0.2 add cross-provider consensus.

3. **Cross-condition DE should use pseudobulk + DESeq2, not single-cell DE.** The single-cell community converged on this in 2023-2024. Naive cell-level Wilcoxon tests overstate significance because of pseudo-replication. Lattice should default to pseudobulk when sample IDs are present in `.obs`.

4. **Shifted log normalization, not CPM/CP10k.** Theis lab `sc-best-practices` 2025 recommendation: use the median count depth, not arbitrary 10^4 or 10^6 targets. Lattice defaults should match.

5. **E2B is mature enough to be the default execution layer.** 15M sandbox sessions/month as of March 2025 means this is no longer a bet — it's the standard.

6. **The Claude Agent SDK is the right framework** for this specific use case, not LangChain/LangGraph. Less abstraction, hooks for safety, model-locked to Claude (which is fine because you want Claude for biology specifically — it leads in the cell-typing benchmarks).

7. **You will be tempted to let the LLM write Scanpy code directly.** Don't. The constrained primitive library is the moat. The Substack issue #40 horror story (LLM-generated code silently reversed case/control, paper went to peer review with wrong conclusions) is the cautionary tale that explains the whole architecture.

---

*End of build spec v0.1*
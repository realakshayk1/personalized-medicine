# PRD: Lattice — The AI Bench Scientist for Single-Cell Biology

**Category:** AI Personalized Medicine (YC RFS) — picks-and-shovels layer
**One-liner:** Lattice is an AI agent that lets bench biologists analyze their own single-cell and bulk RNA-seq data through conversation, without writing Python or R. We are Cursor for biology data analysis.
**Version:** 0.4 — pre-YC application (folds in prereq checklists + companion AGENTS.md)
**Last updated:** May 22, 2026
**Team:** Two undergraduates — UF CS (technical) + Penn LSM (clinical exposure + computational bio)
**YC timeline:** Apply Early Decision (track launched Sept 2025) for Summer 2027 batch after May 2027 graduation
**Status:** No code written yet. This PRD guides the first commit.

**Companion documents:**

- **`BUILD_SPEC.md`** — concrete architecture, 12-week day-by-day plan, primitive registry design, eval methodology. Read this after the PRD if you're about to write code.
- **`AGENTS.md`** — the operating manual for coding agents (Claude Code, Cursor, Codex). Always loaded into agent context. Captures the architectural invariants, code style, and "what to decide vs. ask" rules.
- **`docs/DECISIONS.md`** — ADRs (architecture decision records) for every meaningful tech choice. Add an ADR when you change anything that contradicts the AGENTS.md or BUILD_SPEC.
- **`docs/PRIMITIVES.md`** — the authoritative spec for every primitive in the library.
- **`docs/EVAL.md`** — benchmark methodology, current scores, regression history.

If a human is asking "what are we building and why?", they read this PRD. If an agent is asking "how do I write code for this project?", it reads AGENTS.md. Don't merge them.

---

## TL;DR for a YC partner

> Every paper in personalized medicine, every cancer immunotherapy trial, every CAR-T optimization starts with single-cell RNA-seq. The Benchling 2026 AI Report says 89% of scientists already use AI copilots — but the #1 reason their workflows fail is fragmented data and the gap between bench biologists and the Python/R tools (Scanpy, Seurat) that actually analyze the data. Most biology PhDs don't code well, and the ones who do spend half their time fighting Scanpy instead of doing biology. We're building an agent that ingests an AnnData/h5ad file, lets the user ask questions in English, runs the analysis with auditable code, and returns publication-quality figures with the code, statistics, and figures any reviewer can verify. Two undergrads. The LSM co-founder is the buyer persona. The first 50 customers are at Penn, UF, and the broader academic network. We sell per-seat to labs and per-team to small biotechs. We're not in the regulatory path of clinical products, so we can ship in months, not years.

---

## 0. How this PRD is different from v0.2

I previously recommended an MD-sign-off rare-disease diagnostic agent. That was wrong for this team. Two reasons:

1. **No clinical co-founder, no committed advisor.** That product requires geneticist authority. You don't have it yet. Don't build a product whose moat you can't construct.
2. **Building alongside school.** Clinical-product timelines (BAAs, IRB, sign-off rosters, KOL recruitment) require full-time founders. A B2B research-tools product can be built and sold by undergrads during school.

I also raised concerns about misrepresenting your school status. Don't. **YC's Early Decision track (launched September 2025) is the legitimate path**: apply this fall while still enrolled, get accepted, defer until after May 2027 graduation, do the Summer 2027 batch. Sneha Sivakumar and Anushka Nijhawan (Spur) did exactly this and raised $4.5M. The narrative writes itself: "We applied early, finished school, and now we're full-time." YC partners respect that and the network respects it.

The product below is designed to be the strongest possible thing two undergrads can ship by fall 2026 — not the strongest possible thing in personalized medicine.

---

## 1. The wedge

### 1.1 Why this beats other options for this team

| Option | Fit | Why |
|---|---|---|
| **Clinical AI (diagnostic agent, MTB copilot)** | ❌ | Needs MD credentials/sign-off; HIPAA; long sales cycles. Wrong for undergrads. |
| **Consumer health AI** | ❌ | Crowded (Clara, Nori, Juno, HealthEx, Anthropic itself). Distribution-heavy; needs full-time. |
| **N-of-1 therapeutics** | ❌ | Wet lab, FDA, capital-intensive. |
| **AI agent for bench biology data analysis** | ✅ | LSM co-founder *is* the buyer; Penn + UF lab network is the first 100 users; no regulatory blocker; ships in months; YC has funded analogues (Alchemy P26 for image analysis); $$ market is real (Benchling = $7B+ company built on adjacent pain) |
| **AI for biotech drug discovery** | ⚠️ | Crowded (Variant Bio, Reticular YC, Strand AI YC, Angstrom AI YC, 10x Science YC). PhD-heavy founder norm. Undergrads can play but won't win the funded-by-OrbiMed game. |

The bench-biology data-analysis wedge has a specific structural property: **the user is the buyer, the buyer has no procurement, and the buyer is a current classmate or PI of yours**. That's the only setup where an undergrad team can plausibly hit 100 paying users in 12 months without quitting school.

### 1.2 The specific pain we're solving

The Benchling 2026 AI Report (n=100 biotech/biopharma orgs, Nov 2025 survey):

- 89% of scientists already use AI copilots as their default first stop
- Literature review (76% adoption), protein structure prediction (71%), scientific reporting (66%), target identification (58%) have all broken out of pilot
- **Generative design, biomarker analysis, ADME stay stuck — because the underlying data is scattered, incomplete, and hard to validate**
- "Build what differentiates, buy what scales" — labs are buying tools that handle workflow infrastructure

The specific bottleneck we attack: **a biologist with a single-cell RNA-seq experiment has to either (a) learn enough Python/R + Scanpy/Seurat to do it themselves over weeks, or (b) hand off to a bioinformatician with a 4–8 week queue.** Neither is acceptable. Scanpy and Seurat are powerful but their learning curve makes biologists hate themselves. The new GPU pipelines (RAPIDS-singlecell, ScaleSC) are even more inaccessible.

Two real data points:
- The PLOS Computational Biology piece "Ten Simple Rules for Biologists Learning to Program" — the "how to exit vim" meme literally appears as a figure
- The Substack horror story (issue #40) where an LLM-generated RNA-seq script silently reversed treatment and control groups; the researcher submitted to a journal and only caught it in peer review

The latter is the existential risk to a naïve "ChatGPT analyzes your data" pitch and it shapes our architecture decisively (see §3).

### 1.3 Why now

- **Foundation models can write production-grade analysis code** with the right scaffolding (Claude Code SDK validates this every day)
- **AnnData / h5ad has won as the universal single-cell file format** — one format we can target
- **Scanpy 2.0** (2025) cleaned the Python API; modern Seurat 5 is more programmatic
- **Anthropic launched Claude for Life Sciences in October 2025** with direct Benchling, PubMed, 10x Genomics, and BioRender integrations — the platform layer exists to build on
- **The Nov 2025 Benchling survey is your validation slide** — biopharma is spending on AI right now
- **Alchemy (YC Spring 2026)** is doing the analogous thing for microscopy image analysis. Their existence proves the model. Their existence does NOT cover single-cell or bulk transcriptomics — that's our lane.

---

## 2. The product

### 2.1 What it is

Lattice is a desktop + web app where a bench biologist drags in an `.h5ad`, `.rds`, or raw 10x output, and then has a conversation:

> **User**: "I have a CITE-seq dataset from murine TILs, post-anti-PD1 treatment vs control. Find me T-cell subpopulations that are differentially expanded in responders."
>
> **Lattice**: [generates a plan: QC → normalize → integrate → cluster → annotate using PanglaoDB/CellTypist → DE across conditions → visualize] [executes step-by-step, showing code, log output, and figures] [returns a summary: "Three T-cell subpopulations are significantly expanded in responders: CD8+ effector memory (cluster 4, 3.2× expansion, p=0.0008), exhausted CD8+ (cluster 7, p=0.012), and Treg (cluster 12, p=0.03). The full DE table, UMAP, and reproducible script are below."]
>
> **User**: "Re-run cluster 4 vs cluster 7 with stricter QC and show me the top 20 markers as a heatmap."
>
> **Lattice**: [does it]

The output is always three things together:
1. **A figure or table** suitable for a paper or thesis chapter
2. **The code that produced it** (annotated, reproducible, in the user's choice of Python/Scanpy or R/Seurat)
3. **A verification panel**: the QC numbers, the cell counts at each step, the sanity checks that should make a reviewer trust this

### 2.2 What we are NOT building

- A generic ChatGPT wrapper. The story about LLM-generated code silently reversing case/control groups is the cautionary tale. **We do not let the model write whatever it wants.** Every analysis runs through a constrained set of validated, version-pinned pipelines.
- A new Scanpy or Seurat. Those are tools we orchestrate, not replace.
- A LIMS/ELN. We integrate with Benchling and eLabNext; we don't compete.
- A cloud HPC platform. We run locally on the user's laptop or their lab's existing compute (with optional cloud burst).
- A clinical product. Lattice is Research Use Only. No HIPAA, no FDA.

### 2.3 First wedge: single-cell + bulk RNA-seq for academic labs and small biotechs

Reasons for the narrow start:
- **Highest pain density**: scRNA-seq adoption has exploded; every immunology, cancer, neuroscience, and stem-cell lab now generates this data
- **Standardized inputs**: 10x Genomics + AnnData covers ~80% of new datasets
- **Concrete benchmarks**: we can prove correctness against published datasets where the answer is known (Tabula Muris, Human Cell Atlas, immune-checkpoint trial data)
- **Your network**: every Penn LSM thesis lab and every UF biology PI is a potential design-partner

Expansion path:
1. scRNA-seq → bulk RNA-seq (simpler, larger user base of older labs)
2. Add ATAC-seq, CUT&RUN, spatial transcriptomics (Visium, Xenium, MERFISH)
3. Add CRISPR screen analysis (MAGeCK orchestration)
4. Add flow cytometry as an FCS ingester (compete with OMIQ on the AI-conversation angle)
5. Eventually become "the agent that does anything your bioinformatician would do"

We do NOT do clinical, regulated, or patient-data work in v1, v2, or v3. That door reopens only if/when you add a clinical co-founder.

---

## 3. Architecture

### 3.1 Why constrained tool-use beats free-form code generation

The core insight: a biology graduate student trusts a paper that uses Seurat or Scanpy with the right parameters far more than they trust a model that writes its own pipeline. Reviewers are the same. **Lattice's job is not to write novel bioinformatics code. It's to orchestrate validated tools and explain what they did.**

Concrete design rules:
- The model selects from a **curated, version-pinned library of analysis primitives** (e.g., `scanpy.pp.normalize_total`, `scanpy.tl.leiden`, `scanpy.tl.rank_genes_groups`). The model fills parameters; it does not write the function bodies.
- Every primitive has **input/output schemas** — an AnnData with required `.obs` columns in, an AnnData with named added columns out. Schema validation catches case/control swaps before they happen.
- **All "judgment" calls are explicit and surfaced to the user**: which normalization, which integration method, which clustering resolution. Defaults are evidence-based (citing benchmarks), but never silent.
- **Sanity checks run automatically**: cell counts before/after filtering, library-size distributions, mitochondrial-gene fractions, batch-effect diagnostics. These appear in the verification panel.
- **The reasoning trace is part of the output**, in plain language, so a reviewer or PI can read what was done without reading the code.

### 3.2 Stack

```
┌──────────────────────────────────────────────────────────────┐
│  USER INTERFACE                                                │
│  • Web app (Next.js + TypeScript)                              │
│  • Desktop wrapper (Tauri) for local file access + privacy     │
│  • Chat + plan + figure + code + verification panel layout     │
└──────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  ORCHESTRATOR (Claude Opus 4.7 + Sonnet 4.6)                   │
│  • Planner agent: decomposes user question into a plan        │
│  • Executor agent: selects primitives, fills parameters       │
│  • Verifier agent: runs sanity checks, summarizes results     │
│  Built on Anthropic Agent SDK + Claude Code SDK patterns      │
└──────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  PRIMITIVE LIBRARY (version-pinned, schema-validated)         │
│  • Scanpy 2.x (Python)                                         │
│  • Seurat 5.x (R, via reticulate / Rserve)                    │
│  • RAPIDS-singlecell (GPU, optional)                           │
│  • DESeq2, edgeR (bulk DE)                                     │
│  • MAGeCK (CRISPR)                                             │
│  • CellTypist, PanglaoDB (annotation)                          │
│  • Harmony, scVI, BBKNN (integration)                          │
└──────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  EXECUTION                                                     │
│  • Local: user's machine (Python/R subprocess)                 │
│  • Burst: user's AWS/GCP credentials for big datasets          │
│  • State: AnnData/Seurat objects pickled + versioned locally   │
└──────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  OUTPUTS                                                       │
│  • Publication-grade figures (matplotlib/scanpy/seurat themes) │
│  • Annotated, runnable script (.py or .R)                      │
│  • Verification report (auto-generated markdown)               │
│  • Full provenance: every parameter, every version, every seed │
└──────────────────────────────────────────────────────────────┘
```

### 3.3 The trust layer is the product

This deserves its own section because it is the differentiator vs. "ChatGPT but for biology."

- **Reproducibility**: every analysis exports as a script that runs end-to-end without Lattice. The user owns their work.
- **Provenance**: every figure has a SHA referencing the exact code, parameters, and input data hash.
- **Auditability**: the verification panel surfaces sanity-check failures (mitochondrial-gene fraction abnormally high, batch effect not corrected, cell counts inconsistent across steps) so a PI doing thesis review can catch problems.
- **Citation suggestions**: when Lattice uses Harmony, it surfaces the citation. When it uses CellTypist, it surfaces the model version and reference paper. The biologist's "Methods" section writes itself.
- **Failure modes are loud**: if the model is uncertain which method to use, it asks. If the data violates an assumption (e.g., too few cells per condition for reliable DE), it warns.

### 3.4 What we build in 12 weeks vs. 6 months

**12 weeks (summer 2026 prototype)**:
- Web app + Tauri wrapper
- One workflow well: scRNA-seq from `.h5ad` → QC → cluster → annotate → DE → publication figures
- 10 core primitives validated against ground-truth datasets
- Verification panel
- 20 design-partner users (Penn + UF labs)

**6 months (by EOY 2026)**:
- Bulk RNA-seq workflow
- ATAC-seq workflow
- Multi-dataset integration
- Team / lab account features (shared notebooks)
- 100 paying users; 5 paying labs at $50–200/mo per seat

---

## 4. Users and buyers

| Persona | What they want | What they pay | How we reach |
|---|---|---|---|
| **Grad student / postdoc analyzing scRNA-seq for their thesis** | "Just give me the figure" — bypasses the bioinformatician queue | Personal $20–50/mo, or absorbed by PI | Penn LSM network; UF biology; cold outreach to Twitter/Bluesky bio-academia |
| **PI of a 5–15 person wet lab** | Free her bioinformatician for harder problems; let bench scientists self-serve | Lab license $100–300/seat/mo, 5–10 seats | KOL referrals from design partners |
| **Small biotech (Seed–Series A)** | No dedicated bioinformatics team; runs everything through 1 overworked person | Team license $500–2K/mo | Outbound to YC Biotech batch; reverse-network from your PI advisors |
| **Core facility at a university** | Bottleneck reduction for the queue of incoming projects | Site license $5–25K/yr | After 6 month traction with PIs |
| **Future: pharma R&D team** | Self-serve analysis without going through formal bioinformatics ticket | Enterprise $50K+/yr | Year 2+ after credibility built |

**Primary GTM**: bottom-up viral spread through grad students, with PI/lab licenses as the monetization gate. Mirror Notion, Linear, Cursor.

---

## 5. Competitive landscape (May 2026)

| Competitor | What they do | Where they fall short |
|---|---|---|
| **OMIQ** (BD-owned) | Cloud flow cytometry analysis, cytometry-specific | Cytometry only; not scRNA-seq; not LLM-agent |
| **Trailmaker (10x Genomics)** | Vendor-tied scRNA-seq analysis | 10x customers only; rigid pipeline; not conversational |
| **Galaxy Project** | Open-source workflow runner | UI is brutal; not agent; not AI-native |
| **Latch Bio** | Cloud bioinformatics platform | Engineer-targeted; not biologist-facing |
| **Benchling** | ELN/LIMS + AI copilot integrations | Workflow infra, not analysis; their Anthropic copilot only works if your Benchling is populated |
| **Sphinx Bio** | AI for biotech data | Enterprise pharma sales motion; not bottom-up |
| **Alchemy (YC Spring 2026)** | AI agent for microscopy image analysis | Different modality; adjacent, not direct |
| **Strand AI / 10x Science / Reticular (YC W26)** | Foundation models for biology data | Different layer — they're predicting biology; we're analyzing measurements |
| **Anthropic Claude for Life Sciences (Oct 2025)** | Platform layer | We build on top |
| **Generic ChatGPT / Claude.ai** | Code generation | No primitive validation, no provenance, the case/control-reversal horror story |

**Defensible position**: the only conversational, biologist-targeted, validation-first AI agent for single-cell + bulk transcriptomics. Schema-validated primitive library is the technical moat. Penn + UF lab network is the GTM moat. Two-undergrad team and the timing means we move faster than incumbents.

The real threat is **10x Genomics extending Trailmaker into an AI agent** or **Benchling pushing their Claude integration into analysis territory**. Mitigation: own the bottom-up biology grad-student market they don't reach; integrate with them where they push down; be acquisition-attractive within 2–3 years if we win the user love.

---

## 6. The 12-week plan (summer 2026)

This is the entire plan that the next commit needs to support.

### Weeks 1–2: Foundations
- Set up monorepo (Next.js + FastAPI + Tauri scaffold)
- Stand up Claude Code SDK orchestrator
- Build primitive registry pattern: 5 Scanpy primitives, fully schema-validated
- Ingest `.h5ad` and write to a versioned local store
- One end-to-end pipeline: load → QC → normalize → cluster → return

### Weeks 3–4: First workflow
- 10 Scanpy primitives total (PCA, UMAP, Leiden, rank_genes_groups, etc.)
- Annotation via CellTypist
- Differential expression
- Auto-generated figures (matplotlib + scanpy plotting)
- Verification panel v1

### Weeks 5–6: Design partner alpha
- Onboard 3 Penn LSM labs and 2 UF biology labs
- Run their real datasets end-to-end
- Daily user-feedback loop
- Build the "explain what was done" output

### Weeks 7–8: The polish that wins
- Reproducibility: export script that runs standalone
- Figure styling: publication defaults that don't look AI-generated
- Verification: surface every sanity check
- One-click "compare to published method" for benchmark datasets

### Weeks 9–10: Second workflow + monetization
- Bulk RNA-seq workflow (DESeq2 wrapper)
- Stripe + $20/mo personal tier
- Lab license tier (multi-seat)

### Weeks 11–12: YC application prep
- 20+ active users, 5+ paying
- 3 video testimonials from Penn/UF
- One published paper or thesis that credits Lattice
- Early Decision application submitted (deadline likely late October / early November 2026)

---

## 7. Milestones and metrics

| When | Milestone | Target |
|---|---|---|
| Wk 4 | First end-to-end scRNA-seq workflow | Lattice produces a UMAP + DE table on Tabula Muris that matches the published paper figures |
| Wk 8 | 10 active design-partner users | Each runs Lattice on a real dataset, gives feedback |
| Wk 12 | 25 users, 5 paying | $100–500 MRR; YC Early Decision submitted |
| Mo 6 (EOY 2026) | 100 users, 15 paying labs | $3–8K MRR; published paper with Lattice in Methods |
| Mo 9 (Mar 2027) | YC interview round | Demo: live dataset analyzed in 5 min |
| May 2027 | Graduate, defer YC to Summer 2027 | Both founders full-time |
| Aug 2027 | YC Summer 2027 batch | 200+ users, $20K+ MRR, Series Seed open |

---

## 8. Risks and mitigations

| Risk | Severity | Mitigation |
|---|---|---|
| **LLM-generated code reverses case/control or silently corrupts analysis** | Critical | Constrained primitive library; schema validation; sanity-check verification panel; we never generate analysis code freeform — only orchestrate validated primitives |
| **10x Genomics extends Trailmaker into a real agent** | High | Be acquisition-attractive within 24 months; own bottom-up grad-student love they don't have |
| **Benchling extends their Anthropic copilot into analysis** | Medium | Their wedge is workflow infra (ELN), not analysis. Their copilot is enterprise-sold; ours is grad-student-bottom-up. Different distribution |
| **Anthropic ships a "Claude for Life Sciences" analysis agent themselves** | Medium | Vertical depth beats horizontal: we build the biology-specific primitive library and validation infra Anthropic won't |
| **Two undergrads, both in school, can't ship enough** | High | Use YC Early Decision path; aggressive scope cut to one workflow; both founders commit 25 hr/wk minimum during semester, full-time summer |
| **Biology grad students won't pay** | Medium | They won't personally — PIs will. Lab license at $100–300/seat is in the discretionary spend zone (under PI signing authority) |
| **GPU costs / inference costs make per-seat economics ugly** | Low–Medium | Most work runs locally; LLM cost is the variable; budget $5–15 in inference cost per active user/month; gross margin still 70%+ at $50/seat |
| **One founder loses momentum mid-school** | High | Both equally vested; clear weekly cadence; YC Early Decision creates the external commitment |
| **Existing OSS gets "good enough" with native Claude/GPT integrations** | Medium | The validation layer (sanity checks, primitive schemas, provenance) is real technical IP. The product wins on trust, not raw LLM capability |

---

## 9. Why a YC partner should care

The pitch in one paragraph:

> Every breakthrough in personalized medicine — every CAR-T optimization, every immune-checkpoint biomarker discovery, every rare-disease therapeutic target — starts with someone running scRNA-seq or bulk RNA-seq analysis in Scanpy or Seurat. 89% of those people now reach for an AI copilot first. But generic LLMs hallucinate, generate code that silently reverses case and control groups, and produce results that fail peer review. Lattice is the validated, biologist-targeted, conversation-first agent that closes the gap. We're building Cursor for biology data analysis. Two founders, both technical with deep domain exposure, shipping fast, with the Penn and UF lab network as our first 100 customers. We don't need MDs, BAAs, or FDA approval — we sell software to scientists who buy software, and we are the user persona we're building for.

The wedge isn't "AI personalized medicine" in the most literal sense — but it's the **picks-and-shovels layer that every personalized-medicine company needs**. Tempus, Foundation, Variant Bio, and Mendra all run scRNA-seq workflows. We can sell to all of them eventually. Today we sell to the academic labs that train the next generation of personalized-medicine scientists.

---

## 10. Prerequisites before the first commit

This section is the gate. No code lands in the repo until every box here is checked. Treat it as a contract between the two founders and the coding agents that will help build Lattice.

### 10.1 Customer-discovery validation (Week -2)

| # | Task | Owner | Done when |
|---|---|---|---|
| 1 | Interview 20 biology grad students and postdocs at Penn LSM, UF biology, and Bluesky/Twitter bio-academia | LSM | Notes filed in `docs/research/interviews.md` with verbatim quotes. Specific questions: How often do you run scRNA-seq? How long does it take? Who do you ask for help? Would you pay $29/mo personal? Would your PI pay $99/seat? |
| 2 | Interview 5 PIs about discretionary lab software spend | LSM | Notes filed. What do they currently buy? What do they reject? What does their procurement process look like? |
| 3 | Interview 3 small-biotech (Seed–Series A) computational leads | LSM | Notes filed. Test the $500–2K/mo team-license appetite. |
| 4 | Interview 2 core-facility directors | LSM | Notes filed. Year-2 GTM validation only — do not let this distract Week 1 build. |
| 5 | Run Tabula Muris end-to-end through Claude Opus + Claude Agent SDK manually, no Lattice wrapper | CS | Notebook in `docs/research/day_zero_eval.ipynb`. Document specifically: where does raw Claude fail? Hallucinated function names? Bad parameters? Missing reproducibility? This is the calibration that justifies the primitive registry. |
| 6 | Apply to YC Early Decision (Sept/Oct 2026 deadline) | Both | Application submitted. Even rejection forces the pitch to crystallize. |

**Gate rule:** if fewer than 4 of these 6 produce signal that the wedge is real, do not start the build. Re-stress-test in v0.5 of this PRD instead.

### 10.2 Technical preflight (Week -1)

Concrete environmental and knowledge checks. Both founders must be able to answer yes to every line before the first commit.

| # | Check | Owner |
|---|---|---|
| 1 | Tabula Muris loads in Scanpy on both laptops without error | Both |
| 2 | Both founders have manually clustered the PBMC 3k dataset in Scanpy and can explain every step | Both |
| 3 | Anthropic API key issued, billing set up, ≥$200 development credit available | CS |
| 4 | E2B account created; "claude" template sandbox launched at least once; an `h5ad` uploaded and inspected via Scanpy inside the sandbox | CS |
| 5 | Theis lab `sc-best-practices.org` chapters read: QC, Normalization, HVG, Integration, Annotation, DE | LSM |
| 6 | Repo ownership confirmed: CS owns `apps/desktop`, `apps/orchestrator`, infra; LSM owns `packages/primitives`, `eval/`, biology correctness | Both |
| 7 | Git workflow agreed: feature branches, PRs even for solo work, Conventional Commits, weekly sync cadence given school schedule | Both |
| 8 | Both founders have read `AGENTS.md` and `BUILD_SPEC.md` end-to-end | Both |

### 10.3 Benchmark datasets sourced (Week -1)

The eval suite is part of the prereq, not an afterthought. Without ground truth, there is no product.

| # | Dataset | Source | Used for |
|---|---|---|---|
| 1 | Tabula Muris (mouse cell atlas, 100K cells, multi-tissue) | tabula-muris.ds.czbiohub.org | Easy baseline; multi-tissue annotation |
| 2 | PBMC 3k (10x Genomics canonical tutorial) | 10x Genomics | Floor test; most-replicated dataset in field |
| 3 | Human Cell Atlas immune subset (T-cell subtypes) | data.humancellatlas.org | Fine-grained annotation precision |
| 4–8 | 5 published treatment-vs-control papers from past 12 months (immune checkpoint, CAR-T, drug perturbation) | GEO accessions filed in `eval/benchmarks/README.md` | DE accuracy vs published results |
| 9–10 | 2 "trap" datasets with known QC failure modes (high mito, doublet contamination, batch effect) | Curated from public scRNA-seq with known artifacts | Sanity-check recall (must flag, not silently process) |

Done when all 10 datasets are downloaded, hashed, and stored under `eval/benchmarks/` with metadata in `eval/benchmarks/README.md`.

### 10.4 Sketch deliverables (Week 0)

Before week-1 of the 12-week build, produce:

- **The primitive registry decorator** in Python — implement 3 primitives end-to-end (`calculate_qc_metrics`, `filter_cells_basic`, `normalize_total_log1p`) against a test h5ad. Verify schema validation actually catches a deliberate `condition` column swap. This is the load-bearing experiment for the whole architecture.
- **A planner-prompt prototype** — system prompt + 3 few-shot examples that turn user prompts into ordered primitive invocations. Test against 5 example user queries. Refine until 4/5 produce sensible plans.
- **Frontend wireframes** — chat panel + plan panel + figure panel + code panel layout. Figma or paper. Show to 3 grad-student interviewees from Step 1 above for reaction.

### 10.5 The kill-criterion

If after the preflight any of the following is true, **stop and replan**:

- Customer-discovery interviews reveal grad students/PIs would not pay even $29/mo (price discovery, not feature discovery)
- Day-0 eval shows raw Claude Opus already produces acceptable Scanpy output without scaffolding — the primitive library moat is weaker than expected, rethink positioning
- Either founder discovers they cannot allocate 25 hr/week minimum during the school semester
- A direct competitor (Trailmaker AI, Benchling analysis copilot) ships the same product surface in the meantime

A stop here saves months of wasted work. The plan is in service of the goal, not the other way around.

---

## Appendix A — The honest answer about the school + YC question

You asked whether to misrepresent your school status to YC. The honest answer is no, and here's why it's also the strategically correct answer:

1. **YC Early Decision (launched Sept 2025) is the legitimate path for your situation.** Apply this fall, get the YES, defer to Summer 2027 batch, finish school in May 2027. This is what Sneha Sivakumar and Anushka Nijhawan at Spur did — they applied to YC in fall 2023, finished school May 2024, joined Summer 2024 batch, raised $4.5M.
2. **YC partners ask follow-up logistics questions in interviews.** Who's running point during finals? When do you ship? Where will you be in October? Lying once leads to a cascade you can't sustain over a 30-minute interview.
3. **The YC network has long memory.** A reputation hit at age 22 closes doors at 26 when you'd otherwise be in a strong position.
4. **The story is actually better if you're honest about being undergrads.** "We applied Early Decision, here's the traction we built during school, we're full-time after May 2027" reads as confident and planful. Pretending to defer reads as desperate.

The product PRD above is designed around this honest timeline: ship hard during junior/senior year, hit YC Early Decision deadline in fall 2026, batch in summer 2027.

---

## Appendix B — Why I abandoned the v0.2 rare-disease wedge for this team

A previous version of this PRD recommended a clinical rare-disease diagnostic agent with MD geneticist sign-off. That wedge had real merit but is wrong for two undergrads without a clinical co-founder:

1. **Buyer credibility**: clinical geneticists don't buy software from undergrads. Period.
2. **HIPAA / BAA / SOC 2** infrastructure takes a year and a compliance hire you can't afford.
3. **MD sign-off roster** requires recruiting practicing geneticists at advisory equity — none of them will commit to two undergrads pre-traction.
4. **Sales cycle**: 6–12 months at each children's hospital. Impossible to land while in school.

The bench-biology wedge is **strictly more probable to ship and get traction in your timeline**, and it doesn't foreclose the clinical pivot later (with a clinical co-founder added in years 2–3, Lattice's data infrastructure becomes the foundation for clinical-grade products).

---

## Appendix C — Sources and evidence base

**Market and adoption**:
- Benchling 2026 Biotech AI Report (n=100 biopharma orgs, Nov 2025 survey): 89% of scientists use AI copilots; #1 failure mode is data quality
- 6 ways AI reshaped scientific software in 2025 (R&D World, Dec 2025): Claude for Life Sciences, Benchling/PubMed/10x/BioRender integrations
- "The problem with an AI computational biologist isn't the AI" (Scaling Biotech, Substack, Nov 2024): the verification bottleneck argument
- Lab automation market: $9B in 2025, $24B projected 2035 (Grand View Research / SNS Insider)

**Technical foundations**:
- Scanpy benchmark (bioRxiv, Oct 2025): 1.3M-cell datasets, Seurat takes 3 hr vs RAPIDS 20 min — scale problem is real
- scDown (PMC, 2025) and SeuratExtend (GigaScience, July 2025): proliferation of "make this easier" tools shows the demand
- SCAPE (bioRxiv, Nov 2025): AI-driven scRNA-seq platforms are emerging — confirms direction
- ScaleSC (Bioinformatics Advances, 2025): GPU pipelines exist; biologists still can't use them

**The case for biologist-targeted UX**:
- "Ten Simple Rules for Biologists Learning to Program" (PLOS Comp Bio, 2017): the canonical "how to exit vim" pain
- "Data challenges of biomedical researchers" (PMC, n=129): 47% of researchers blame lack of programming training as the data-analysis blocker
- Issue #40 from Sequence and Destroy (Substack): the LLM RNA-seq case/control reversal — why naïve LLM code generation is dangerous

**Competitive context**:
- Alchemy (YC P26): "automates image analysis for life science researchers" — direct proof of model for adjacent modality
- Strand AI, 10x Science, Reticular (YC W26): foundation models for biology
- OMIQ vs FlowJo: cloud-based, AI-augmented analysis tools work — they just don't cover transcriptomics
- Trailmaker (10x Genomics): vendor-tied; biggest direct threat if they pivot to agent

**YC pathway**:
- TechCrunch (Sept 2025): YC Early Decision launch announcement
- Spur (Summer 2024 batch): Sivakumar & Nijhawan, applied as students, deferred, raised $4.5M
- YC blog (Jared Friedman): how to start biotech from a university context

---

*End of PRD v0.4*
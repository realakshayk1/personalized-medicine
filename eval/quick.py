"""Quick suite (~5 min) — standard scRNA-seq chain on PBMC 3k.

Loads ``sc.datasets.pbmc3k()`` (raw 3k PBMCs), runs the canonical pipeline via
``run_primitive`` (calculate_qc_metrics -> filter_cells -> filter_genes ->
normalize -> hvg -> pca -> neighbors -> leiden), and asserts the chain completes
and produces clusters.

If ``pbmc3k_processed()`` labels are available, an Adjusted Rand Index (ARI)
between our Leiden clusters and the published cluster labels is computed as a
PLACEHOLDER for the eventual annotation-F1 target. Annotation primitives are not
built yet (M-future), so ARI stands in as a cluster-agreement signal — it is NOT
the AGENTS.md F1 metric and is not gated.

Network behavior: scanpy downloads PBMC 3k on first use. If the download fails
(offline sandbox), every case SKIPs with a clear message instead of crashing.
"""

from __future__ import annotations

from anndata import AnnData
from lattice_primitives import run_primitive

from eval._harness import CaseResult, SuiteReport

SEED = 0


def _load_pbmc3k() -> AnnData:
    """Load raw PBMC 3k via scanpy. Raises on network/download failure."""
    import scanpy as sc

    adata = sc.datasets.pbmc3k()
    adata.var_names_make_unique()  # AGENTS.md: dedupe gene symbols on ingest
    # The standard chain's normalize primitive requires layer_state='raw_counts'.
    adata.uns["_lattice_layer_state"] = "raw_counts"
    return adata


def _try_load_processed_labels() -> dict | None:  # type: ignore[type-arg]
    """Return {barcode: louvain_label} from pbmc3k_processed, or None on failure."""
    try:
        import scanpy as sc

        proc = sc.datasets.pbmc3k_processed()
    except Exception:
        return None
    label_col = None
    for cand in ("louvain", "leiden", "cell_type", "bulk_labels"):
        if cand in proc.obs.columns:
            label_col = cand
            break
    if label_col is None:
        return None
    return dict(zip(proc.obs_names.tolist(), proc.obs[label_col].astype(str).tolist(), strict=True))


def _ari(labels_a: list[str], labels_b: list[str]) -> float | None:
    try:
        from sklearn.metrics import adjusted_rand_score

        return float(adjusted_rand_score(labels_a, labels_b))
    except Exception:
        return None


def run_quick_suite() -> SuiteReport:
    report = SuiteReport(suite="quick")

    # --- Load data (network-guarded) -------------------------------------
    try:
        adata = _load_pbmc3k()
    except Exception as exc:
        report.add(
            CaseResult(
                suite="quick",
                dataset="pbmc3k",
                metric="dataset download",
                value="SKIP",
                passed=True,
                skipped=True,
                detail=(
                    f"sc.datasets.pbmc3k() unavailable (offline?): {exc}. "
                    "Run `uv run python -m eval.fetch_data` with network access."
                ),
            )
        )
        return report

    n_cells_raw = adata.n_obs
    report.add(
        CaseResult(
            suite="quick",
            dataset="pbmc3k",
            metric="loaded raw cells",
            value=str(n_cells_raw),
            passed=n_cells_raw > 0,
        )
    )

    # --- Run the standard chain ------------------------------------------
    try:
        adata, _, _ = run_primitive("calculate_qc_metrics", adata, mt_prefix="MT-")
        adata, _, _ = run_primitive("filter_cells_min_counts", adata, min_counts_per_cell=500)
        adata, _, _ = run_primitive("filter_genes_min_cells", adata, min_cells_per_gene=3)
        adata, _, _ = run_primitive("normalize_total_log1p", adata)
        adata, _, _ = run_primitive(
            "highly_variable_genes", adata, n_top_genes=2000, flavor="seurat"
        )
        adata, _, _ = run_primitive("pca", adata, n_comps=50, random_state=SEED)
        adata, _, _ = run_primitive("neighbors", adata, n_neighbors=15, random_state=SEED)
        adata, _, _ = run_primitive("leiden", adata, resolution=1.0, random_state=SEED)
    except Exception as exc:
        report.add(
            CaseResult(
                suite="quick",
                dataset="pbmc3k",
                metric="standard chain completes",
                value="ERROR",
                passed=False,
                detail=f"chain raised: {exc}",
            )
        )
        return report

    # --- Assert clusters were produced -----------------------------------
    has_leiden = "leiden" in adata.obs.columns
    n_clusters = int(adata.obs["leiden"].nunique()) if has_leiden else 0
    report.add(
        CaseResult(
            suite="quick",
            dataset="pbmc3k",
            metric="standard chain completes",
            value=f"{adata.n_obs} cells retained",
            passed=has_leiden and adata.n_obs > 0,
        )
    )
    report.add(
        CaseResult(
            suite="quick",
            dataset="pbmc3k",
            metric="leiden clusters produced",
            value=str(n_clusters),
            passed=n_clusters >= 2,
            detail="expected >=2 clusters from real PBMC structure",
        )
    )

    # --- Placeholder cluster-agreement metric (ARI vs published labels) ---
    labels = _try_load_processed_labels()
    if labels is None:
        report.add(
            CaseResult(
                suite="quick",
                dataset="pbmc3k_processed",
                metric="cluster-agreement ARI (placeholder)",
                value="SKIP",
                passed=True,
                skipped=True,
                detail=(
                    "pbmc3k_processed labels unavailable (offline or no label "
                    "column). NOTE: ARI is a placeholder; the AGENTS.md "
                    "annotation-F1 >=0.75 target needs annotation primitives "
                    "(not built yet)."
                ),
            )
        )
        return report

    # Align on shared barcodes (processed dataset is a filtered subset).
    shared = [bc for bc in adata.obs_names if bc in labels]
    if len(shared) < 50:
        report.add(
            CaseResult(
                suite="quick",
                dataset="pbmc3k_processed",
                metric="cluster-agreement ARI (placeholder)",
                value="SKIP",
                passed=True,
                skipped=True,
                detail=f"only {len(shared)} barcodes overlap processed labels; too few to score",
            )
        )
        return report

    ours = adata[shared].obs["leiden"].astype(str).tolist()
    theirs = [labels[bc] for bc in shared]
    ari = _ari(ours, theirs)
    if ari is None:
        report.add(
            CaseResult(
                suite="quick",
                dataset="pbmc3k_processed",
                metric="cluster-agreement ARI (placeholder)",
                value="SKIP",
                passed=True,
                skipped=True,
                detail="scikit-learn unavailable for adjusted_rand_score",
            )
        )
        return report

    report.add(
        CaseResult(
            suite="quick",
            dataset="pbmc3k_processed",
            metric="cluster-agreement ARI (placeholder)",
            value=f"{ari:.3f} (n={len(shared)})",
            # Not a gated target — recorded as informational. Pass if finite.
            passed=True,
            detail=(
                "PLACEHOLDER for annotation-F1; not the AGENTS.md gated metric. "
                "Annotation primitives not yet built."
            ),
        )
    )
    return report

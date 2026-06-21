"""Trap suite — injected failure modes that MUST be caught (100% recall gate).

Each trap builds a synthetic AnnData whose data is engineered to trigger a
specific primitive sanity check, runs the primitive via ``run_primitive``, and
asserts that a warning of the expected severity fired with an expected substring
in its message.

A trap PASSES only when a matching warning is found. The suite as a whole passes
only when EVERY trap passes — AGENTS.md mandates 100% trap-detection recall.

Traps implemented:
    high_mito          calculate_qc_metrics  -> 'warn'  (mitochondrial content)
    not_raw_counts     normalize_total_log1p -> 'error' (max>15 / not raw counts)
    under_clustering   leiden                -> 'warn'  (single/degenerate cluster)
    hvg_overrequest    highly_variable_genes -> 'warn'  (n_top_genes > n_genes)
    overaggressive_filter  filter_cells_min_counts -> 'error' (>80% removed)
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import scanpy as sc
import scipy.sparse
from anndata import AnnData
from lattice_primitives import run_primitive

from eval._harness import CaseResult, SuiteReport

# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _base_counts(n_obs: int, n_vars: int, seed: int, lam: float = 10.0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.poisson(lam=lam, size=(n_obs, n_vars)).astype(np.float32)


def _make_adata(
    counts: np.ndarray,
    *,
    n_mt: int = 10,
    layer_state: str | None = None,
) -> AnnData:
    """Wrap a count matrix in an AnnData with MT-prefixed genes + condition."""
    n_obs, n_vars = counts.shape
    var_names = [f"MT-GENE{i:03d}" for i in range(n_mt)] + [
        f"GENE{i:03d}" for i in range(n_mt, n_vars)
    ]
    obs_names = [f"CELL_{i:04d}" for i in range(n_obs)]
    condition = ["treated"] * (n_obs // 2) + ["control"] * (n_obs - n_obs // 2)

    adata = AnnData(X=scipy.sparse.csr_matrix(counts))
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obs["condition"] = condition
    if layer_state is not None:
        adata.uns["_lattice_layer_state"] = layer_state
    return adata


def _make_high_mito() -> AnnData:
    """90% of counts come from the 10 MT genes -> median pct_counts_mt >> 20%."""
    n_obs, n_vars = 200, 100
    rng = np.random.default_rng(7)
    counts = np.zeros((n_obs, n_vars), dtype=np.float32)
    counts[:, :10] = rng.poisson(lam=500, size=(n_obs, 10)).astype(np.float32)
    counts[:, 10:] = rng.poisson(lam=2, size=(n_obs, n_vars - 10)).astype(np.float32)
    return _make_adata(counts)


def _make_not_raw_counts() -> AnnData:
    """Counts engineered so log1p(normalized) max exceeds 15 -> 'error' fires.

    The check is ``max(adata.X) > 15`` after normalize_total (target_sum=None,
    i.e. scale each cell to the MEDIAN total) + log1p. To exceed 15, a single
    gene must hold close to the full per-cell total after rescaling, and that
    total must be very large: log1p(x) > 15 requires x > ~3.27e6.

    So we give every cell a single dominant gene holding ~5e6 counts plus a
    sparse background. The median total is ~5e6; after scaling to it each cell's
    dominant gene stays ~5e6, and log1p(5e6) ~ 15.4 > 15. This is plainly not
    raw single-cell counts, which is exactly what the 'error' guard protects.
    """
    n_obs, n_vars = 120, 80
    rng = np.random.default_rng(11)
    counts = rng.poisson(lam=1, size=(n_obs, n_vars)).astype(np.float64)
    # One dominant gene per cell carrying ~5e6 counts.
    dominant = rng.integers(low=10, high=n_vars, size=n_obs)
    counts[np.arange(n_obs), dominant] = 5_000_000.0
    return _make_adata(counts.astype(np.float32), layer_state="raw_counts")


def _make_structureless_log_normalized() -> AnnData:
    """Uniform-noise data with no separable structure -> leiden under-clusters.

    Built through the real preprocessing chain so leiden's requires are met, but
    the underlying signal is structureless, so leiden collapses to one (or a
    single dominant) cluster at the default resolution.
    """
    n_obs, n_vars = 300, 400
    rng = np.random.default_rng(3)
    # Identical Poisson rate for every gene/cell: no population structure.
    counts = rng.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)
    adata = _make_adata(counts, n_mt=0, layer_state="raw_counts")
    # Normalize + log1p directly (mirrors normalize_total_log1p output state).
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)
    adata.uns["_lattice_layer_state"] = "log_normalized"
    return adata


def _make_hvg_overrequest() -> AnnData:
    """Small gene count so n_top_genes (5000) far exceeds available genes."""
    n_obs, n_vars = 200, 60
    counts = _base_counts(n_obs, n_vars, seed=5, lam=8.0)
    adata = _make_adata(counts, n_mt=0, layer_state="raw_counts")
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)
    adata.uns["_lattice_layer_state"] = "log_normalized"
    return adata


def _make_overaggressive_filter() -> AnnData:
    """Most cells have tiny counts so a high min_counts removes >80% of cells."""
    n_obs, n_vars = 200, 50
    rng = np.random.default_rng(13)
    # Most cells ~1-2 counts/gene (total well under threshold), a few high.
    counts = rng.poisson(lam=1, size=(n_obs, n_vars)).astype(np.float32)
    counts[:10, :] = rng.poisson(lam=200, size=(10, n_vars)).astype(np.float32)
    adata = _make_adata(counts, n_mt=0)
    # filter_cells_min_counts requires obs['total_counts'] from QC metrics.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    return adata


# ---------------------------------------------------------------------------
# Trap specification
# ---------------------------------------------------------------------------


class Trap:
    """One trap: a fixture + a primitive call + the expected warning."""

    def __init__(
        self,
        name: str,
        primitive: str,
        builder: Callable[[], AnnData],
        params: dict,  # type: ignore[type-arg]
        expect_severity: str,
        expect_substring: str,
    ) -> None:
        self.name = name
        self.primitive = primitive
        self.builder = builder
        self.params = params
        self.expect_severity = expect_severity
        self.expect_substring = expect_substring.lower()


TRAPS: list[Trap] = [
    Trap(
        name="high_mito",
        primitive="calculate_qc_metrics",
        builder=_make_high_mito,
        params={"mt_prefix": "MT-"},
        expect_severity="warn",
        expect_substring="mitochondrial",
    ),
    Trap(
        name="not_raw_counts",
        primitive="normalize_total_log1p",
        builder=_make_not_raw_counts,
        params={},
        expect_severity="error",
        expect_substring="not raw integer counts",
    ),
    Trap(
        name="under_clustering",
        primitive="leiden",
        builder=_make_structureless_log_normalized,
        params={"resolution": 0.1, "random_state": 0},
        expect_severity="warn",
        # Either "only one cluster" or "largest cluster holds .. %" qualifies;
        # both contain "cluster".
        expect_substring="cluster",
    ),
    Trap(
        name="hvg_overrequest",
        primitive="highly_variable_genes",
        builder=_make_hvg_overrequest,
        params={"n_top_genes": 5000, "flavor": "seurat"},
        expect_severity="warn",
        expect_substring="only",
    ),
    Trap(
        name="overaggressive_filter",
        primitive="filter_cells_min_counts",
        builder=_make_overaggressive_filter,
        params={"min_counts_per_cell": 5000},
        expect_severity="error",
        expect_substring="removed",
    ),
]


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def _run_trap(trap: Trap) -> CaseResult:
    """Execute one trap; for the under_clustering trap, neighbors must run first."""
    adata = trap.builder()

    # The leiden trap needs the full dim-reduce prefix so leiden's requires pass.
    if trap.primitive == "leiden":
        adata, _, _ = run_primitive(
            "highly_variable_genes", adata, n_top_genes=200, flavor="seurat"
        )
        adata, _, _ = run_primitive("pca", adata, n_comps=20, random_state=0)
        adata, _, _ = run_primitive("neighbors", adata, n_neighbors=15, random_state=0)

    try:
        _, warnings, _ = run_primitive(trap.primitive, adata, **trap.params)
    except Exception as exc:  # a raised exception is NOT a caught sanity warning
        return CaseResult(
            suite="trap",
            dataset=trap.name,
            metric=f"{trap.primitive} {trap.expect_severity}",
            value="ERROR",
            passed=False,
            detail=f"{trap.primitive} raised instead of warning: {exc}",
        )

    matched = [
        w
        for w in warnings
        if w.severity == trap.expect_severity and trap.expect_substring in w.message.lower()
    ]
    passed = len(matched) > 0
    if passed:
        value = f"{trap.expect_severity} fired"
        detail = matched[0].message
    else:
        value = "NOT CAUGHT"
        fired = ", ".join(f"{w.severity}:{w.message[:40]}" for w in warnings) or "none"
        detail = (
            f"expected severity='{trap.expect_severity}' containing "
            f"'{trap.expect_substring}'; got: {fired}"
        )

    return CaseResult(
        suite="trap",
        dataset=trap.name,
        metric=f"{trap.primitive} {trap.expect_severity}",
        value=value,
        passed=passed,
        detail=detail,
    )


def run_trap_suite() -> SuiteReport:
    """Run all traps. The caller enforces the 100%-recall gate via all_passed."""
    report = SuiteReport(suite="trap")
    for trap in TRAPS:
        report.add(_run_trap(trap))
    return report

"""Primitive: volcano.

Renders a volcano plot for one group from a rank_genes_groups result. Scanpy has
no built-in volcano plot, so this is drawn with raw matplotlib from the marker
result: x = log2 fold change, y = -log10(adjusted p-value) for the chosen group.

Figure transport convention (see plot_umap.py): plot primitives append a base64
PNG data-URI to adata.uns["_lattice_figures"]; run_primitive lifts it into
PrimitiveResult.figures and removes it from .uns.

Server-context rules (AGENTS.md): matplotlib Agg backend is forced, plt.show()
is never called, and figures are always closed.

NOTE (biology-semantics, provisional — see docs/BIOLOGY_DEFAULTS.md):
    The significance thresholds (lfc_threshold=1.0 log2FC, pval_threshold=0.05
    adjusted p) used only for point coloring are provisional and FLAGGED FOR
    LSM-FOUNDER REVIEW. They are display thresholds, not statistical cutoffs, and
    do not alter the underlying rank_genes_groups result.
"""

from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")  # headless backend — must be set before pyplot import

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from anndata import AnnData  # noqa: E402
from loguru import logger  # noqa: E402

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_  # noqa: E402
from lattice_primitives.registry import primitive, sanity_check  # noqa: E402

FIGURES_UNS_KEY = "_lattice_figures"

# ---------------------------------------------------------------------------
# Primitive registration
# ---------------------------------------------------------------------------

_REQUIRES = AnnDataState(
    layer_state="log_normalized",
    # rank_genes_groups (uns) is guarded in the body.
)
_PRODUCES = AnnDataState(
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "group": ParamSpec.model_validate(
        {
            "type": "str",
            "default": None,
            "description": (
                "The group (cluster label) to plot. Defaults to the first group "
                "in the rank_genes_groups result."
            ),
        }
    ),
    "lfc_threshold": ParamSpec.model_validate(
        {
            "type": "float",
            "default": 1.0,
            "description": (
                "Display-only |log2 fold change| threshold above which points are "
                "highlighted. Does not alter the underlying test."
            ),
            "min_value": 0.0,
        }
    ),
    "pval_threshold": ParamSpec.model_validate(
        {
            "type": "float",
            "default": 0.05,
            "description": (
                "Display-only adjusted-p threshold below which points are "
                "highlighted. Does not alter the underlying test."
            ),
            "min_value": 0.0,
            "max_value": 1.0,
        }
    ),
}


@primitive(
    name="volcano",
    category="plot",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_rank_genes_groups",
    row_modifying=False,
)
def volcano(
    adata: AnnData,
    group: str | None = None,
    lfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
) -> AnnData:
    """Render a volcano plot for one group and stash it as a PNG.

    Parameters
    ----------
    adata:
        AnnData with .uns['rank_genes_groups'] (run 'rank_genes_groups' first).
    group:
        Group label to plot; defaults to the first group in the result.
    lfc_threshold:
        Display-only |log2FC| highlight threshold.
    pval_threshold:
        Display-only adjusted-p highlight threshold.

    Returns
    -------
    The same AnnData, with a base64 PNG data-URI appended to
    adata.uns["_lattice_figures"].
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "volcano requires marker results (adata.uns['rank_genes_groups']). "
            "Run the 'rank_genes_groups' primitive before 'volcano'."
        )

    rgg = adata.uns["rank_genes_groups"]
    available_groups = list(rgg["names"].dtype.names)
    if group is None:
        group = available_groups[0]
    elif group not in available_groups:
        raise ValueError(
            f"volcano: group '{group}' not found in rank_genes_groups result. "
            f"Available groups: {available_groups}"
        )

    lfc = np.asarray(rgg["logfoldchanges"][group], dtype=float)
    pvals_adj = np.asarray(rgg["pvals_adj"][group], dtype=float)

    # -log10(p); clip p at a tiny floor so p=0 does not produce inf.
    p_floor = np.nextafter(0, 1)
    neg_log10_p = -np.log10(np.clip(pvals_adj, p_floor, 1.0))

    finite = np.isfinite(lfc) & np.isfinite(neg_log10_p)
    significant = finite & (np.abs(lfc) >= lfc_threshold) & (pvals_adj < pval_threshold)

    fig, ax = plt.subplots(figsize=(6, 5))
    try:
        ax.scatter(
            lfc[finite & ~significant],
            neg_log10_p[finite & ~significant],
            s=8,
            c="lightgray",
            label="ns",
        )
        ax.scatter(
            lfc[significant],
            neg_log10_p[significant],
            s=10,
            c="crimson",
            label="significant",
        )
        ax.axhline(-np.log10(pval_threshold), color="gray", ls="--", lw=0.7)
        ax.axvline(lfc_threshold, color="gray", ls="--", lw=0.7)
        ax.axvline(-lfc_threshold, color="gray", ls="--", lw=0.7)
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10(adjusted p-value)")
        ax.set_title(f"Volcano — group '{group}' (one-vs-rest)")
        ax.legend(loc="upper right", frameon=False, fontsize=8)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        data_uri = "data:image/png;base64," + base64.b64encode(buf.read()).decode("ascii")
    finally:
        plt.close(fig)

    adata.uns.setdefault(FIGURES_UNS_KEY, []).append(data_uri)
    logger.info(
        f"volcano: rendered group='{group}' "
        f"({int(significant.sum())} highlighted, {len(data_uri)} bytes b64)"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("volcano")
def _check_volcano(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate that a figure was produced."""
    warnings: list[Warning_] = []

    figures = adata_after.uns.get(FIGURES_UNS_KEY, [])
    if not figures:
        warnings.append(
            Warning_(
                severity="error",
                message="volcano did not produce a figure.",
                primitive="volcano",
            )
        )

    return warnings

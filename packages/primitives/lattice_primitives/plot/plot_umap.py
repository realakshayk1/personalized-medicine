"""Primitive: plot_umap.

Renders the UMAP embedding colored by a categorical .obs column (default the
Leiden clusters) and emits it as a PNG.

Figure transport convention:
    Plot primitives do not add columns to the AnnData. Instead they append a
    base64-encoded PNG data-URI to adata.uns["_lattice_figures"]. The registry's
    run_primitive helper lifts this list into PrimitiveResult.figures and removes
    it from .uns (so figures do not bloat the persisted AnnData). The orchestrator
    streams figures to the frontend FigurePanel via the execute SSE channel.

Server-context rules (AGENTS.md): the matplotlib Agg backend is forced and the
figure is always closed to avoid leaking memory. plt.show() is never called.
"""

from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")  # headless backend — must be set before pyplot import

import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402
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
)
_PRODUCES = AnnDataState(
    # Emits a figure (via uns["_lattice_figures"], lifted by run_primitive) but
    # adds no schema-tracked AnnData state.
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "color": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "leiden",
            "description": (
                "The .obs column used to color points. Defaults to the Leiden cluster "
                "labels. Must be an existing .obs column."
            ),
        }
    ),
}


@primitive(
    name="plot_umap",
    category="plot",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_calculate_qc_metrics",
    row_modifying=False,
)
def plot_umap(
    adata: AnnData,
    color: str = "leiden",
) -> AnnData:
    """Render the UMAP embedding colored by `color` and stash it as a PNG.

    Parameters
    ----------
    adata:
        AnnData with adata.obsm['X_umap'] present (run 'umap' first).
    color:
        .obs column to color by (default the Leiden clusters).

    Returns
    -------
    The same AnnData, with a base64 PNG data-URI appended to
    adata.uns["_lattice_figures"].
    """
    if "X_umap" not in adata.obsm:
        raise ValueError(
            "plot_umap requires a UMAP embedding (adata.obsm['X_umap']). "
            "Run the 'umap' primitive before 'plot_umap'."
        )
    if color not in adata.obs.columns:
        raise ValueError(
            f"plot_umap: color column '{color}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    fig, ax = plt.subplots(figsize=(6, 5))
    try:
        sc.pl.umap(adata, color=color, ax=ax, show=False)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        data_uri = "data:image/png;base64," + base64.b64encode(buf.read()).decode("ascii")
    finally:
        plt.close(fig)

    adata.uns.setdefault(FIGURES_UNS_KEY, []).append(data_uri)
    logger.info(f"plot_umap: rendered UMAP colored by '{color}' ({len(data_uri)} bytes b64)")
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("plot_umap")
def _check_plot_umap(
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
                message="plot_umap did not produce a figure.",
                primitive="plot_umap",
            )
        )

    return warnings

"""Primitive: heatmap.

Renders a Scanpy expression heatmap of the top marker genes per group and emits
it as a PNG. Uses sc.pl.rank_genes_groups_heatmap, which reads the marker result
in adata.uns['rank_genes_groups'] directly.

Figure transport convention (see plot_umap.py): plot primitives append a base64
PNG data-URI to adata.uns["_lattice_figures"]; run_primitive lifts it into
PrimitiveResult.figures and removes it from .uns.

Server-context rules (AGENTS.md): matplotlib Agg backend is forced, plt.show()
is never called, and figures are always closed.
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
    # rank_genes_groups (uns) and the groupby column are guarded in the body.
)
_PRODUCES = AnnDataState(
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "groupby": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "leiden",
            "description": (
                "The .obs column whose groups label the heatmap. Should match the "
                "groupby used for rank_genes_groups."
            ),
        }
    ),
    "n_genes": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 10,
            "description": "Number of top marker genes per group to display.",
            "min_value": 1,
        }
    ),
}


@primitive(
    name="heatmap",
    category="plot",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_rank_genes_groups",
    row_modifying=False,
)
def heatmap(
    adata: AnnData,
    groupby: str = "leiden",
    n_genes: int = 10,
) -> AnnData:
    """Render a marker-gene heatmap and stash it as a PNG.

    Parameters
    ----------
    adata:
        AnnData with .uns['rank_genes_groups'] (run 'rank_genes_groups' first)
        and the `groupby` column in .obs.
    groupby:
        .obs column defining heatmap row groups (default the Leiden clusters).
    n_genes:
        Top markers per group to show.

    Returns
    -------
    The same AnnData, with a base64 PNG data-URI appended to
    adata.uns["_lattice_figures"].
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "heatmap requires marker results (adata.uns['rank_genes_groups']). "
            "Run the 'rank_genes_groups' primitive before 'heatmap'."
        )
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"heatmap: groupby column '{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # show=False returns the axes dict; we grab the active figure to save it.
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_genes,
        groupby=groupby,
        show=False,
        show_gene_labels=True,
    )
    fig = plt.gcf()
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        data_uri = "data:image/png;base64," + base64.b64encode(buf.read()).decode("ascii")
    finally:
        plt.close("all")

    adata.uns.setdefault(FIGURES_UNS_KEY, []).append(data_uri)
    logger.info(
        f"heatmap: rendered top {n_genes} markers/group across "
        f"groupby='{groupby}' ({len(data_uri)} bytes b64)"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("heatmap")
def _check_heatmap(
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
                message="heatmap did not produce a figure.",
                primitive="heatmap",
            )
        )

    return warnings

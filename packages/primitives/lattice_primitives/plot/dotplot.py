"""Primitive: dotplot.

Renders a Scanpy dotplot of the top marker genes per group and emits it as a
PNG. Marker genes are pulled from adata.uns['rank_genes_groups'] (top `n_genes`
per group) so the dotplot summarizes the differential-expression result.

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
                "The .obs column whose groups label the dotplot rows. Should "
                "match the groupby used for rank_genes_groups."
            ),
        }
    ),
    "n_genes": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 5,
            "description": "Number of top marker genes per group to display.",
            "min_value": 1,
        }
    ),
}


def _top_markers_per_group(adata: AnnData, n_genes: int) -> list[str]:
    """Collect the top `n_genes` marker gene names per group, de-duplicated."""
    names = adata.uns["rank_genes_groups"]["names"]
    groups = list(names.dtype.names)
    ordered: list[str] = []
    seen: set[str] = set()
    for g in groups:
        for gene in list(names[g])[:n_genes]:
            gene_str = str(gene)
            if gene_str not in seen:
                seen.add(gene_str)
                ordered.append(gene_str)
    return ordered


@primitive(
    name="dotplot",
    category="plot",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_rank_genes_groups",
    row_modifying=False,
)
def dotplot(
    adata: AnnData,
    groupby: str = "leiden",
    n_genes: int = 5,
) -> AnnData:
    """Render a dotplot of top markers per group and stash it as a PNG.

    Parameters
    ----------
    adata:
        AnnData with .uns['rank_genes_groups'] (run 'rank_genes_groups' first)
        and the `groupby` column in .obs.
    groupby:
        .obs column defining dotplot rows (default the Leiden clusters).
    n_genes:
        Top markers per group to show.

    Returns
    -------
    The same AnnData, with a base64 PNG data-URI appended to
    adata.uns["_lattice_figures"].
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "dotplot requires marker results (adata.uns['rank_genes_groups']). "
            "Run the 'rank_genes_groups' primitive before 'dotplot'."
        )
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"dotplot: groupby column '{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    var_names = _top_markers_per_group(adata, n_genes)
    if not var_names:
        raise ValueError("dotplot: no marker genes found in adata.uns['rank_genes_groups'].")

    # return_fig=True yields a DotPlot object; .savefig writes its composite figure.
    dp = sc.pl.dotplot(
        adata,
        var_names=var_names,
        groupby=groupby,
        return_fig=True,
        show=False,
    )
    try:
        buf = io.BytesIO()
        dp.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        data_uri = "data:image/png;base64," + base64.b64encode(buf.read()).decode("ascii")
    finally:
        plt.close("all")

    adata.uns.setdefault(FIGURES_UNS_KEY, []).append(data_uri)
    logger.info(
        f"dotplot: rendered {len(var_names)} markers across "
        f"groupby='{groupby}' ({len(data_uri)} bytes b64)"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("dotplot")
def _check_dotplot(
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
                message="dotplot did not produce a figure.",
                primitive="dotplot",
            )
        )

    return warnings

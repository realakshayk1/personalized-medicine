"""Primitives: filter_cells_min_counts + filter_genes_min_cells.

This module registers TWO primitives that are run sequentially as a coupled
pair in the standard scRNA-seq pipeline. They are split (rather than combined)
for cleaner provenance logging — each step has its own input/output hash.

Design note — why two primitives instead of one:
    sc.pp.filter_cells and sc.pp.filter_genes operate on different axes (obs
    vs var) and have independent thresholds. Keeping them separate allows the
    planner to use only one (e.g., apply cell filtering without gene filtering)
    and makes the provenance log show exactly how many cells and how many genes
    were removed in each step. AGENTS.md §"How to add a primitive" permits a
    "tightly-coupled pair" in one file; we split for maximum auditability.

Requires:
    obs column 'total_counts' (from calculate_qc_metrics).

Params:
    filter_cells_min_counts: min_counts_per_cell: int = 500
    filter_genes_min_cells:  min_cells_per_gene: int = 3
"""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_
from lattice_primitives.registry import primitive, sanity_check

# ---------------------------------------------------------------------------
# filter_cells_min_counts
# ---------------------------------------------------------------------------

_REQ_CELLS = AnnDataState(
    obs_columns=["total_counts"],
)
_PROD_CELLS = AnnDataState(
    immutable_obs_columns=["condition"],
)
_PARAMS_CELLS = {
    "min_counts_per_cell": ParamSpec.model_validate({
        "type": "int",
        "default": 500,
        "description": (
            "Minimum total UMI count for a cell to be retained. "
            "Cells with fewer counts are removed as low-quality. "
            "Default 500 follows the Theis lab sc-best-practices guidance."
        ),
        "min_value": 1,
    }),
}


@primitive(
    name="filter_cells_min_counts",
    category="preprocess",
    requires=_REQ_CELLS,
    produces=_PROD_CELLS,
    params=_PARAMS_CELLS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_filter_cells",
    row_modifying=True,  # removes cells — obs count changes
)
def filter_cells_min_counts(
    adata: AnnData,
    min_counts_per_cell: int = 500,
) -> AnnData:
    """Filter cells with fewer than `min_counts_per_cell` total UMI counts.

    Wraps sc.pp.filter_cells. This is row-modifying: cells that do not meet
    the threshold are removed from the AnnData object.

    Parameters
    ----------
    adata:
        AnnData post calculate_qc_metrics (requires obs['total_counts']).
    min_counts_per_cell:
        Minimum count threshold. Cells below this are dropped.

    Returns
    -------
    Filtered AnnData with low-count cells removed.
    """
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    n_after = adata.n_obs
    logger.info(
        f"filter_cells_min_counts: removed {n_before - n_after} cells "
        f"({n_before} -> {n_after}, threshold={min_counts_per_cell})"
    )
    return adata


@sanity_check("filter_cells_min_counts")
def _check_filter_cells(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    warnings: list[Warning_] = []
    # Use _n_obs_before injected by run_primitive to handle in-place mutations
    # where adata_before IS adata_after (same object after scanpy in-place filter).
    n_before = params.get("_n_obs_before", adata_before.n_obs)
    n_after = adata_after.n_obs
    pct_removed = 100.0 * (n_before - n_after) / max(n_before, 1)

    if pct_removed > 80.0:
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    f"{pct_removed:.1f}% of cells removed (>{80}%) with "
                    f"min_counts_per_cell={params.get('min_counts_per_cell', 500)}. "
                    "This suggests the threshold is too aggressive. "
                    "Verify the input data and consider lowering the threshold."
                ),
                primitive="filter_cells_min_counts",
            )
        )
    elif pct_removed > 50.0:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"{pct_removed:.1f}% of cells removed (>50%) with "
                    f"min_counts_per_cell={params.get('min_counts_per_cell', 500)}. "
                    "Verify this reflects expected low-quality cell fraction."
                ),
                primitive="filter_cells_min_counts",
            )
        )

    return warnings


# ---------------------------------------------------------------------------
# filter_genes_min_cells
# ---------------------------------------------------------------------------

_REQ_GENES = AnnDataState()  # no required obs columns — genes don't need QC metrics
_PROD_GENES = AnnDataState(
    immutable_obs_columns=["condition"],
)
_PARAMS_GENES = {
    "min_cells_per_gene": ParamSpec.model_validate({
        "type": "int",
        "default": 3,
        "description": (
            "Minimum number of cells in which a gene must be detected (count > 0) "
            "for the gene to be retained. Default 3 follows Theis lab guidance."
        ),
        "min_value": 1,
    }),
}


@primitive(
    name="filter_genes_min_cells",
    category="preprocess",
    requires=_REQ_GENES,
    produces=_PROD_GENES,
    params=_PARAMS_GENES,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_filter_cells",
    row_modifying=False,  # removes genes (var), not cells (obs)
)
def filter_genes_min_cells(
    adata: AnnData,
    min_cells_per_gene: int = 3,
) -> AnnData:
    """Filter genes detected in fewer than `min_cells_per_gene` cells.

    Wraps sc.pp.filter_genes. This removes rarely-expressed genes that
    contribute noise without biological signal. Obs count is unchanged;
    only adata.n_vars decreases.

    Parameters
    ----------
    adata:
        AnnData object.
    min_cells_per_gene:
        Minimum cell count for a gene to be retained.

    Returns
    -------
    AnnData with low-detection genes removed.
    """
    n_vars_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    n_vars_after = adata.n_vars
    logger.info(
        f"filter_genes_min_cells: removed {n_vars_before - n_vars_after} genes "
        f"({n_vars_before} -> {n_vars_after}, threshold={min_cells_per_gene} cells)"
    )
    return adata


@sanity_check("filter_genes_min_cells")
def _check_filter_genes(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    warnings: list[Warning_] = []
    n_before = params.get("_n_vars_before", adata_before.n_vars)
    n_after = adata_after.n_vars
    pct_removed = 100.0 * (n_before - n_after) / max(n_before, 1)

    if pct_removed > 80.0:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"{pct_removed:.1f}% of genes removed with "
                    f"min_cells_per_gene={params.get('min_cells_per_gene', 3)}. "
                    "Verify this matches expected sparsity for your dataset type."
                ),
                primitive="filter_genes_min_cells",
            )
        )

    return warnings

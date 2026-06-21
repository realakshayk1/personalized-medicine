"""Primitive: calculate_qc_metrics.

Wraps sc.pp.calculate_qc_metrics to add standard quality-control columns to
.obs and .var. This is the canonical first step in the Theis lab best-practices
pipeline before any filtering or normalization.

Columns added to .obs:
    total_counts         — total UMI/read count per cell
    n_genes_by_counts    — number of genes detected per cell (count > 0)
    pct_counts_mt        — fraction of counts from mitochondrial genes (%)

Columns added to .var:
    mean_counts          — mean count across cells
    n_cells_by_counts    — number of cells in which the gene is detected

Mitochondrial gene detection: genes whose name starts with `mt_prefix`
(default "MT-" for human). Mouse datasets typically use lowercase "mt-".
The primitive auto-detects the lowercase variant if no genes match the
provided prefix.
"""

from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_
from lattice_primitives.registry import primitive, sanity_check

# ---------------------------------------------------------------------------
# Primitive registration
# ---------------------------------------------------------------------------

_REQUIRES = AnnDataState()  # no prerequisites; first primitive in pipeline
_PRODUCES = AnnDataState(
    obs_columns=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    var_columns=["mean_counts", "n_cells_by_counts"],
    # layer_state is passthrough — whatever was in before stays
    immutable_obs_columns=["condition"],  # never let QC clobber condition labels
)
_PARAMS = {
    "mt_prefix": ParamSpec.model_validate({
        "type": "str",
        "default": "MT-",
        "description": "Prefix for mitochondrial gene names. 'MT-' for human; 'mt-' for mouse. Auto-detects lowercase if no genes match.",
    }),
}


@primitive(
    name="calculate_qc_metrics",
    category="qc",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_calculate_qc_metrics",
    row_modifying=False,
)
def calculate_qc_metrics(
    adata: AnnData,
    mt_prefix: str = "MT-",
) -> AnnData:
    """Compute per-cell and per-gene QC metrics via sc.pp.calculate_qc_metrics.

    Parameters
    ----------
    adata:
        AnnData object. `.X` may be sparse or dense.
    mt_prefix:
        Prefix for mitochondrial gene names. Defaults to "MT-" (human).
        The function auto-detects lowercase "mt-" as a fallback for mouse data
        when no genes match the provided prefix.

    Returns
    -------
    AnnData with QC columns added in-place.
    """
    # Determine mitochondrial gene prefix
    prefix = mt_prefix
    mt_mask = adata.var_names.str.startswith(prefix)
    n_mt = mt_mask.sum()

    if n_mt == 0:
        # Auto-detect lowercase for mouse datasets
        lower_prefix = prefix.lower()
        if lower_prefix != prefix:
            mt_mask_lower = adata.var_names.str.startswith(lower_prefix)
            n_mt_lower = mt_mask_lower.sum()
            if n_mt_lower > 0:
                logger.info(
                    f"No genes matched mt_prefix='{prefix}'; "
                    f"auto-detected {n_mt_lower} genes with lowercase '{lower_prefix}'."
                )
                prefix = lower_prefix
                mt_mask = mt_mask_lower

    if mt_mask.sum() == 0:
        logger.warning(
            f"No mitochondrial genes found with prefix '{prefix}'. "
            "pct_counts_mt will be 0.0 for all cells."
        )

    adata.var["mt"] = mt_mask

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # Rename scanpy's "pct_counts_mt" convention is already correct; verify
    # scanpy adds: total_counts, n_genes_by_counts, pct_counts_mt (to obs)
    #              mean_counts, n_cells_by_counts (to var)

    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("calculate_qc_metrics")
def _check_qc_metrics(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    warnings: list[Warning_] = []

    median_pct_mt = float(np.median(adata_after.obs["pct_counts_mt"]))

    if median_pct_mt > 20.0:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Median mitochondrial content is {median_pct_mt:.1f}% (>20%). "
                    "This may indicate low-quality tissue, stressed cells, or poor dissociation. "
                    "Consider stricter mitochondrial filtering."
                ),
                primitive="calculate_qc_metrics",
            )
        )

    return warnings

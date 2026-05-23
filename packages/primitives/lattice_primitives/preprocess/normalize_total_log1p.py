"""Primitive: normalize_total_log1p.

Wraps sc.pp.normalize_total followed by sc.pp.log1p as a documented tightly-
coupled pair per AGENTS.md. These two operations are inseparable in the Theis
lab best-practices pipeline — normalize then immediately log-transform so that
downstream HVG selection, PCA, and clustering operate on a stabilized scale.

Reference: Luecken & Theis, Molecular Systems Biology, 2019
    (citation_key: theis_normalization_2019)

Layer state:
    requires: "raw_counts"  (set by the ingest primitive in the full pipeline)
    produces: "log_normalized"

Modifications to AnnData:
    adata.layers["counts_normalized"] — pre-normalization copy of raw counts
    adata.X                           — replaced with log1p(normalized counts)
    adata.uns["_lattice_layer_state"] — set to "log_normalized"

The `condition` obs column is declared immutable: any primitive that swaps
case/control labels will be detected and rejected by the registry before the
corrupted data can propagate downstream.
"""

from __future__ import annotations

import numpy as np
import scanpy as sc
import scipy.sparse
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_
from lattice_primitives.registry import primitive, sanity_check

# ---------------------------------------------------------------------------
# Primitive registration
# ---------------------------------------------------------------------------

_REQUIRES = AnnDataState(
    layer_state="raw_counts",
)
_PRODUCES = AnnDataState(
    layer_state="log_normalized",
    adds_layers=["counts_normalized"],
    immutable_obs_columns=["condition"],  # case/control protection
)
_PARAMS = {
    "target_sum": ParamSpec.model_validate({
        "type": "float",
        "default": None,
        "description": (
            "Per-cell scaling target for normalize_total. "
            "None (default) uses the median total_counts across cells — "
            "the Theis lab recommendation (Luecken & Theis 2019). "
            "Set to 1e4 for CP10k or 1e6 for CPM."
        ),
    }),
    "exclude_highly_expressed": ParamSpec.model_validate({
        "type": "bool",
        "default": False,
        "description": (
            "If True, excludes the top-1% most highly expressed genes "
            "from the normalization denominator. Recommended for datasets "
            "dominated by a single cell type (e.g., erythrocyte-rich tissue)."
        ),
    }),
}


@primitive(
    name="normalize_total_log1p",
    category="preprocess",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="theis_normalization_2019",
    row_modifying=False,
)
def normalize_total_log1p(
    adata: AnnData,
    target_sum: float | None = None,
    exclude_highly_expressed: bool = False,
) -> AnnData:
    """Normalize counts per cell then apply log1p transformation.

    Step 1 — Save raw counts:
        adata.layers["counts_normalized"] = adata.X.copy()
        (named "counts_normalized" because it will hold the normalized counts
        after the next step; the pre-norm copy is what the orchestrator uses
        for rollback and provenance)

    Step 2 — Normalize per cell:
        sc.pp.normalize_total(adata, target_sum=target_sum,
                              exclude_highly_expressed=exclude_highly_expressed)
        Each cell's total count is scaled to target_sum (or the median).

    Step 3 — Log-transform:
        sc.pp.log1p(adata)
        Applies log(x + 1) element-wise to adata.X.

    Parameters
    ----------
    adata:
        AnnData in raw_counts layer_state. .X must contain raw integer counts.
    target_sum:
        Normalization target. None -> median of raw total counts (recommended).
    exclude_highly_expressed:
        Exclude top-1% genes from normalization denominator.

    Returns
    -------
    AnnData with adata.X = log1p(normalized) and raw copy in layers["counts_normalized"].
    """
    # Save pre-normalization raw counts
    adata.layers["counts_normalized"] = adata.X.copy()

    # Normalize per cell (in-place)
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
        inplace=True,
    )

    # Log1p transform (in-place)
    sc.pp.log1p(adata)

    logger.info(
        f"normalize_total_log1p: target_sum={target_sum}, "
        f"exclude_highly_expressed={exclude_highly_expressed}"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("normalize_total_log1p")
def _check_normalization(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate normalization output.

    Checks:
    1. Max value > 15 after log1p => input may not have been raw counts.
    2. Median normalized total deviates from log1p(target_sum) by > 1.0.
    """
    warnings: list[Warning_] = []

    X_after = adata_after.X
    if scipy.sparse.issparse(X_after):
        max_val = float(X_after.max())
        row_sums = np.asarray(X_after.sum(axis=1)).ravel()
    else:
        X_arr = np.asarray(X_after)
        max_val = float(X_arr.max())
        row_sums = X_arr.sum(axis=1)

    # Check 1: max value
    if max_val > 15.0:
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    f"Maximum value after log1p is {max_val:.2f} (>15). "
                    "This strongly suggests the input was not raw integer counts. "
                    "Verify adata.X contains raw counts before running this primitive."
                ),
                primitive="normalize_total_log1p",
            )
        )

    # Check 2: median total deviation
    target_sum = params.get("target_sum")
    if target_sum is None:
        # Median of original total counts
        if "total_counts" in adata_before.obs.columns:
            raw_median = float(np.median(adata_before.obs["total_counts"]))
        else:
            # Compute from the raw layer snapshot
            X_before = adata_before.X
            if scipy.sparse.issparse(X_before):
                raw_median = float(np.median(np.asarray(X_before.sum(axis=1)).ravel()))
            else:
                raw_median = float(np.median(np.asarray(X_before).sum(axis=1)))
        expected_log = float(np.log1p(raw_median))
    else:
        expected_log = float(np.log1p(target_sum))

    median_total = float(np.median(row_sums))
    deviation = abs(median_total - expected_log)

    if deviation > 1.0:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Median normalized total per cell ({median_total:.2f}) deviates "
                    f"from expected log1p(target_sum)={expected_log:.2f} by {deviation:.2f} "
                    "(threshold: 1.0). Check for extreme outlier cells or incorrect target_sum."
                ),
                primitive="normalize_total_log1p",
            )
        )

    return warnings

"""Primitive: pseudobulk.

Wraps decoupler's ``dc.pp.pseudobulk`` to aggregate single-cell raw counts into
pseudobulk profiles, one row per (sample x group). Pseudobulk is the prerequisite
for cross-condition differential expression with bulk RNA-seq methods (DESeq2):
aggregating cells back to the biological-replicate level restores the sample-level
variance structure those methods assume, avoiding the inflated false positives
that come from treating individual cells as replicates (Squair et al. 2021).

This primitive is STRUCTURALLY UNUSUAL: it AGGREGATES cells, so it returns a NEW
AnnData whose n_obs is (n_samples x n_groups), not the input cell count. Its
.obs_names become the pseudobulk sample IDs ("{sample}_{group}"). It is therefore
declared ``row_modifying=True``.

decoupler API note (decoupler 2.x): the current entrypoint is
``dc.pp.pseudobulk(adata, sample_col, groups_col, layer=None, raw=False,
mode='sum', ...)``. It sums RAW counts per (sample, group) and returns a fresh
AnnData. (The pre-2.0 ``dc.get_pseudobulk`` no longer exists.)

Requires raw counts. ``sample_col`` and ``groups_col`` are PARAMETERS, so they
cannot be declared in the static ``requires.obs_columns``; the body guards on
their presence explicitly and raises ValueError if missing.

Additions / changes to AnnData (a NEW AnnData is returned):
    .X                            — summed raw counts per (sample x group)
    .obs                          — carries forward sample-constant metadata
                                    (e.g. 'condition'), plus decoupler's
                                    'psbulk_cells' / 'psbulk_counts' bookkeeping
    .uns["_lattice_layer_state"]  — set to "pseudobulk_counts"

NOTE (biology-semantics, provisional — FLAGGED FOR LSM-FOUNDER REVIEW):
    - mode='sum' (sum raw counts) is the DESeq2-compatible default. 'mean' is
      NOT appropriate for count-based DE. Provisional default; do not change
      without LSM sign-off.
    - min_cells=10: pseudobulk samples aggregated from fewer than 10 cells are
      dropped as statistically unreliable. The 10-cell floor is a provisional
      default (the field uses anywhere from 5 to 30); flagged for LSM review.
"""

from __future__ import annotations

import decoupler as dc
import numpy as np
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
    # sample_col / groups_col are parameters, not statically known, so they are
    # guarded in the body rather than declared in requires.obs_columns.
)
_PRODUCES = AnnDataState(
    layer_state="pseudobulk_counts",
    # NOTE: pseudobulk REINDEXES rows by sample ID, so the immutable-column guard
    # for "condition" effectively no-ops here (no surviving cell barcodes map
    # across the before/after obs_names). We keep it declared for consistency and
    # documentation; the real case/control protection for this track lives in the
    # explicit, oriented contrast in pseudobulk_deseq2.
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "sample_col": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "sample",
            "description": (
                "The .obs column identifying the biological replicate / sample "
                "each cell came from. Pseudobulk profiles are summed within each "
                "(sample, group). Must be an existing .obs column."
            ),
        }
    ),
    "groups_col": ParamSpec.model_validate(
        {
            "type": "str",
            "default": None,
            "description": (
                "The .obs column defining the cell grouping (e.g. 'leiden' or a "
                "cell-type column), so DE is computed per cell population. If "
                "None, all cells of a sample are aggregated into one profile. "
                "Must be an existing .obs column when provided."
            ),
        }
    ),
    "mode": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "sum",
            "description": (
                "Aggregation mode. 'sum' (default) sums raw counts and is the "
                "only mode appropriate for downstream count-based DESeq2."
            ),
            "choices": ["sum", "mean"],
        }
    ),
    "min_cells": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 10,
            "description": (
                "Minimum number of cells per pseudobulk sample. Profiles "
                "aggregated from fewer cells are dropped. Provisional default."
            ),
            "min_value": 0,
        }
    ),
}


@primitive(
    name="pseudobulk",
    category="preprocess",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",  # scanpy unused here; pin tracked for consistency
    citation_key="badia_decoupler_2022",
    row_modifying=True,
)
def pseudobulk(
    adata: AnnData,
    sample_col: str = "sample",
    groups_col: str | None = None,
    mode: str = "sum",
    min_cells: int = 10,
) -> AnnData:
    """Aggregate single-cell raw counts to pseudobulk profiles via decoupler.

    Parameters
    ----------
    adata:
        AnnData in raw_counts layer_state. .X must contain raw integer counts.
    sample_col:
        .obs column identifying the biological replicate / sample.
    groups_col:
        .obs column defining the cell grouping (e.g. cell type / cluster). None
        aggregates all of a sample's cells into a single profile.
    mode:
        Aggregation mode, "sum" (recommended for DESeq2) or "mean".
    min_cells:
        Drop pseudobulk samples built from fewer than this many cells.

    Returns
    -------
    A NEW AnnData of shape (n_samples x n_groups) with summed raw counts in .X,
    sample-level metadata carried into .obs, and layer_state="pseudobulk_counts".
    """
    if sample_col not in adata.obs.columns:
        raise ValueError(
            f"pseudobulk requires sample_col '{sample_col}' in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}."
        )
    if groups_col is not None and groups_col not in adata.obs.columns:
        raise ValueError(
            f"pseudobulk requires groups_col '{groups_col}' in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}."
        )

    pb = dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=groups_col,
        mode=mode,
    )

    # Filter low-cell pseudobulk samples. decoupler records the per-profile cell
    # count in .obs["psbulk_cells"]; fall back to no filtering if absent.
    n_before = pb.n_obs
    if min_cells > 0 and "psbulk_cells" in pb.obs.columns:
        keep = pb.obs["psbulk_cells"].to_numpy() >= min_cells
        pb = pb[keep].copy()
    n_after = pb.n_obs

    # Carry over the session layer_state tag (decoupler returns a fresh AnnData
    # with empty .uns); the decorator also sets this from produces.layer_state,
    # but we set it explicitly so the value is correct even outside run_primitive.
    pb.uns["_lattice_layer_state"] = "pseudobulk_counts"

    logger.info(
        f"pseudobulk: sample_col='{sample_col}', groups_col='{groups_col}', "
        f"mode='{mode}', min_cells={min_cells} -> {n_after} pseudobulk samples "
        f"({n_before - n_after} dropped) x {pb.n_vars} genes"
    )
    return pb


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("pseudobulk")
def _check_pseudobulk(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate pseudobulk output.

    Checks:
    1. Too few pseudobulk samples for a meaningful DE comparison (< 4).
    2. .X is not integer-valued (DESeq2 requires raw integer counts).
    3. A condition column, if present, lacks >= 2 levels with >= 2 replicates each.
    """
    warnings: list[Warning_] = []

    # 1: replicate count
    if adata_after.n_obs < 4:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Only {adata_after.n_obs} pseudobulk sample(s) remain after "
                    "aggregation/filtering. DESeq2 needs several biological "
                    "replicates per group; results from <4 samples are unreliable."
                ),
                primitive="pseudobulk",
            )
        )

    # 2: integer counts
    X = adata_after.X
    if scipy.sparse.issparse(X):
        x_arr = np.asarray(X.todense())
    else:
        x_arr = np.asarray(X)
    if x_arr.size > 0 and not np.allclose(x_arr, np.round(x_arr)):
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "Pseudobulk .X is not integer-valued. DESeq2 requires summed "
                    "raw integer counts; verify the input was raw counts and "
                    "mode='sum' (not 'mean')."
                ),
                primitive="pseudobulk",
            )
        )

    # 3: condition replication
    if "condition" in adata_after.obs.columns:
        counts = adata_after.obs["condition"].value_counts()
        levels_with_reps = int((counts >= 2).sum())
        if levels_with_reps < 2:
            warnings.append(
                Warning_(
                    severity="warn",
                    message=(
                        "Fewer than two condition levels have >= 2 pseudobulk "
                        "replicates. Cross-condition DE needs replication on both "
                        "sides of the contrast (e.g. >= 2 treated and >= 2 control)."
                    ),
                    primitive="pseudobulk",
                )
            )

    return warnings

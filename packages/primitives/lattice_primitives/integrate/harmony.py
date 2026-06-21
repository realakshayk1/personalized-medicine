"""Primitive: harmony.

Wraps sc.external.pp.harmony_integrate to correct batch effects in the PCA
embedding via Harmony (Korsunsky et al. 2019). Harmony iteratively projects
cells into a shared embedding where batch-specific variation is regressed out
while biological structure is preserved.

Requires pca to have run first (reads adata.obsm['X_pca']). The requires-schema
cannot express an obsm prerequisite, so the body guards on X_pca explicitly. It
also guards that the chosen batch_key column exists in adata.obs.

Additions to AnnData:
    adata.obsm["X_pca_harmony"]  — n_obs x n_pcs batch-corrected embedding

Downstream note: after harmony, the 'neighbors' primitive should be run with
use_rep="X_pca_harmony" instead of the default "X_pca" so the integrated
embedding feeds clustering/UMAP. Changing neighbors to thread this is a future
planner concern and intentionally out of scope here.

Layer state is passthrough (log_normalized in, log_normalized out).
random_state is threaded for reproducibility (AGENTS.md: non-negotiable).
harmony_integrate forwards **kwargs to harmonypy.run_harmony, which accepts
random_state (verified via inspect) — we pass it through directly.

NOTE (biology-semantics, provisional — flag for LSM-founder review):
    batch_key default "batch" is provisional; the correct column is dataset-
    specific (sample, donor, 10x lane, etc.). All harmonypy tuning parameters
    (theta, lamb, sigma) are left at their defaults; whether to expose theta
    (diversity-clustering penalty) as a tunable param is an open question for
    the LSM founder.
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

_REQUIRES = AnnDataState(
    layer_state="log_normalized",
)
_PRODUCES = AnnDataState(
    adds_obsm=["X_pca_harmony"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "batch_key": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "batch",
            "description": (
                "The .obs column identifying batches/samples to integrate over "
                "(e.g. sample, donor, or 10x lane). Cells sharing a value are "
                "treated as one batch."
            ),
        }
    ),
    "random_state": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 0,
            "description": "Random seed. Threaded from the session seed for reproducibility.",
        }
    ),
}


@primitive(
    name="harmony",
    category="integrate",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="korsunsky_harmony_2019",
    row_modifying=False,
)
def harmony(
    adata: AnnData,
    batch_key: str = "batch",
    random_state: int = 0,
) -> AnnData:
    """Batch-correct the PCA embedding via sc.external.pp.harmony_integrate.

    Parameters
    ----------
    adata:
        AnnData in log_normalized layer_state with adata.obsm['X_pca'] present
        (run 'pca' first).
    batch_key:
        Name of the adata.obs column identifying batches to integrate over.
    random_state:
        Seed forwarded to harmonypy.run_harmony for reproducibility.

    Returns
    -------
    AnnData with .obsm['X_pca_harmony'] added in-place.
    """
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "harmony requires a PCA embedding (adata.obsm['X_pca']). "
            "Run the 'pca' primitive before 'harmony'."
        )
    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"harmony requires batch_key='{batch_key}' to be a column in adata.obs, "
            f"but it is missing. Available obs columns: {list(adata.obs.columns)}."
        )

    # harmony_integrate forwards **kwargs to harmonypy.run_harmony, which accepts
    # random_state (verified via inspect.signature). Thread it for reproducibility.
    sc.external.pp.harmony_integrate(
        adata,
        key=batch_key,
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        random_state=random_state,
    )
    n_batches = int(adata.obs[batch_key].nunique())
    logger.info(
        f"harmony: integrated over batch_key='{batch_key}' "
        f"({n_batches} batches, random_state={random_state}) -> X_pca_harmony"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("harmony")
def _check_harmony(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate Harmony integration output.

    Checks:
    1. X_pca_harmony present (error if missing).
    2. X_pca_harmony is all-finite (error on NaN/Inf — degenerate integration).
    3. X_pca_harmony has the same shape as X_pca (error otherwise).
    4. Only one batch level present -> integration is a no-op (warn).
    """
    warnings: list[Warning_] = []

    x_harmony = adata_after.obsm.get("X_pca_harmony")
    if x_harmony is None:
        warnings.append(
            Warning_(
                severity="error",
                message="X_pca_harmony was not written to adata.obsm. Harmony did not run.",
                primitive="harmony",
            )
        )
        return warnings

    x_harmony_arr = np.asarray(x_harmony)
    if not np.all(np.isfinite(x_harmony_arr)):
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "X_pca_harmony contains non-finite values (NaN/Inf). The "
                    "integration may have diverged; check the input PCA and batch labels."
                ),
                primitive="harmony",
            )
        )

    x_pca = adata_after.obsm.get("X_pca")
    if x_pca is not None and np.asarray(x_pca).shape != x_harmony_arr.shape:
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    f"X_pca_harmony shape {x_harmony_arr.shape} does not match X_pca "
                    f"shape {np.asarray(x_pca).shape}. Harmony should preserve dimensionality."
                ),
                primitive="harmony",
            )
        )

    batch_key = str(params.get("batch_key", "batch"))
    if batch_key in adata_after.obs.columns:
        n_batches = int(adata_after.obs[batch_key].nunique())
        if n_batches <= 1:
            warnings.append(
                Warning_(
                    severity="warn",
                    message=(
                        f"Only one batch level present in '{batch_key}'. Harmony "
                        "integration is a no-op when there is nothing to integrate over; "
                        "verify batch_key points at the correct column."
                    ),
                    primitive="harmony",
                )
            )

    return warnings

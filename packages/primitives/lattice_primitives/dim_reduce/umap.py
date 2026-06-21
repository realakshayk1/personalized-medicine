"""Primitive: umap.

Wraps sc.tl.umap to compute a 2D UMAP embedding for visualization. UMAP reads the
neighbor graph built by the 'neighbors' primitive; it does not build its own.

Requires neighbors to have run first (reads adata.uns['neighbors']). The body
guards on this explicitly since the requires-schema cannot express a uns key.

Additions to AnnData:
    adata.obsm["X_umap"]    — n_obs x 2 UMAP coordinates
    adata.uns["umap"]       — UMAP parameters

Layer state is passthrough. random_state is threaded for reproducibility.

Reproducibility note: UMAP via numba can be non-bit-identical across machines even
with a fixed seed (multithreading). Bit-identical eval may require single-threaded
execution. This does not affect within-machine determinism.
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
    adds_obsm=["X_umap"],
    adds_uns=["umap"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "min_dist": ParamSpec.model_validate(
        {
            "type": "float",
            "default": 0.5,
            "description": (
                "Minimum distance between embedded points. Lower values pack clusters "
                "more tightly. 0.5 is the scanpy default."
            ),
            "min_value": 0.0,
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
    name="umap",
    category="dim_reduce",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_umap_neighbors",
    row_modifying=False,
)
def umap(
    adata: AnnData,
    min_dist: float = 0.5,
    random_state: int = 0,
) -> AnnData:
    """Compute a 2D UMAP embedding via sc.tl.umap.

    Parameters
    ----------
    adata:
        AnnData with a neighbor graph (run 'neighbors' first).
    min_dist:
        Minimum distance between embedded points.
    random_state:
        Seed for reproducibility.

    Returns
    -------
    AnnData with .obsm['X_umap'] and .uns['umap'] added in-place.
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "umap requires a neighbor graph (adata.uns['neighbors']). "
            "Run the 'neighbors' primitive before 'umap'."
        )

    sc.tl.umap(
        adata,
        min_dist=min_dist,
        random_state=random_state,
    )
    logger.info(f"umap: min_dist={min_dist}, random_state={random_state}")
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("umap")
def _check_umap(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate UMAP output.

    Checks:
    1. X_umap missing or contains non-finite values.
    """
    warnings: list[Warning_] = []

    x_umap = adata_after.obsm.get("X_umap")
    if x_umap is None:
        warnings.append(
            Warning_(
                severity="error",
                message="UMAP coordinates (X_umap) were not produced.",
                primitive="umap",
            )
        )
    elif not np.all(np.isfinite(np.asarray(x_umap))):
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "X_umap contains non-finite values (NaN/Inf). The neighbor "
                    "graph may be disconnected or degenerate."
                ),
                primitive="umap",
            )
        )

    return warnings

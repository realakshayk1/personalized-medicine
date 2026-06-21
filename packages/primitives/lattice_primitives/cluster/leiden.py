"""Primitive: leiden.

Wraps sc.tl.leiden to partition the neighbor graph into communities (clusters).
Leiden is the field-standard graph-based clustering for scRNA-seq.

Requires neighbors to have run first (reads adata.uns['neighbors']). The body
guards on this explicitly.

Version note (scanpy 1.11): scanpy emits a FutureWarning unless flavor is set,
because the default backend is moving from 'leidenalg' to 'igraph'. We pass
flavor="igraph", n_iterations=2, directed=False — the documented future defaults,
faster and warning-free. Requires the 'igraph' package.

NOTE (biology-semantics, provisional — see docs/BIOLOGY_DEFAULTS.md):
    resolution=1.0 is provisional. Cluster count is highly sensitive to resolution.
    BUILD_SPEC §1.3 step 8 calls for a multi-resolution sweep (0.2/0.5/0.8/1.2)
    that the user picks from; that is a deliberate follow-up.

Columns added to .obs:
    leiden    — categorical cluster labels (key set by key_added, default "leiden")

Layer state is passthrough. random_state is threaded for reproducibility.
"""

from __future__ import annotations

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
    obs_columns=["leiden"],
    adds_uns=["leiden"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "resolution": ParamSpec.model_validate(
        {
            "type": "float",
            "default": 1.0,
            "description": (
                "Clustering resolution. Higher values yield more, smaller clusters. "
                "1.0 is the scanpy default; 0.5 gives cleaner top-level clusters."
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
    name="leiden",
    category="cluster",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="traag_leiden_2019",
    row_modifying=False,
)
def leiden(
    adata: AnnData,
    resolution: float = 1.0,
    random_state: int = 0,
) -> AnnData:
    """Cluster cells via sc.tl.leiden (igraph backend).

    Parameters
    ----------
    adata:
        AnnData with a neighbor graph (run 'neighbors' first).
    resolution:
        Clustering resolution.
    random_state:
        Seed for reproducibility.

    Returns
    -------
    AnnData with categorical .obs['leiden'] cluster labels added in-place.
    """
    if "neighbors" not in adata.uns:
        raise ValueError(
            "leiden requires a neighbor graph (adata.uns['neighbors']). "
            "Run the 'neighbors' primitive before 'leiden'."
        )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=random_state,
        key_added="leiden",
    )
    n_clusters = adata.obs["leiden"].nunique()
    logger.info(
        f"leiden: resolution={resolution}, random_state={random_state} -> {n_clusters} clusters"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("leiden")
def _check_leiden(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate clustering output.

    Checks:
    1. Only one cluster found (under-clustering / resolution too low).
    2. A single cluster contains >95% of cells (effectively no structure).
    3. More clusters than ~ n_obs/10 (over-clustering / resolution too high).
    """
    warnings: list[Warning_] = []

    labels = adata_after.obs["leiden"]
    n_clusters = int(labels.nunique())
    n_obs = adata_after.n_obs

    if n_clusters <= 1:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    "Leiden found only one cluster. The resolution may be too low or "
                    "the data may lack separable structure. Try a higher resolution."
                ),
                primitive="leiden",
            )
        )
    else:
        largest_frac = float(labels.value_counts(normalize=True).iloc[0])
        if largest_frac > 0.95:
            warnings.append(
                Warning_(
                    severity="warn",
                    message=(
                        f"The largest cluster holds {largest_frac * 100:.0f}% of cells. "
                        "Clustering may be under-resolved; consider a higher resolution."
                    ),
                    primitive="leiden",
                )
            )

    if n_clusters > max(10, n_obs // 10):
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Leiden produced {n_clusters} clusters for {n_obs} cells, which "
                    "may indicate over-clustering. Consider a lower resolution."
                ),
                primitive="leiden",
            )
        )

    return warnings

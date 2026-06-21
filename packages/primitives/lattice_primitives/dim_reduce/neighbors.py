"""Primitive: neighbors.

Wraps sc.pp.neighbors to build the k-nearest-neighbor graph on the PCA embedding.
The neighbor graph is the shared substrate for both UMAP and Leiden clustering.

Requires pca to have run first (reads adata.obsm['X_pca']). The requires-schema
cannot express an obsm prerequisite, so the body guards on X_pca explicitly.

Additions to AnnData:
    adata.uns["neighbors"]          — graph parameters
    adata.obsp["distances"]         — sparse kNN distances (documented; not schema-validated)
    adata.obsp["connectivities"]    — sparse weighted graph (documented; not schema-validated)

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
    adds_uns=["neighbors"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "n_neighbors": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 15,
            "description": (
                "Size of the local neighborhood used for graph construction. 15 is "
                "the scanpy default; larger values give more global structure."
            ),
            "min_value": 2,
        }
    ),
    "n_pcs": ParamSpec.model_validate(
        {
            "type": "int",
            "default": None,
            "description": (
                "Number of principal components to use from X_pca. None uses all "
                "available components."
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
    name="neighbors",
    category="dim_reduce",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_umap_neighbors",
    row_modifying=False,
)
def neighbors(
    adata: AnnData,
    n_neighbors: int = 15,
    n_pcs: int | None = None,
    random_state: int = 0,
) -> AnnData:
    """Build the kNN graph via sc.pp.neighbors on the PCA embedding.

    Parameters
    ----------
    adata:
        AnnData with adata.obsm['X_pca'] present (run pca first).
    n_neighbors:
        Local neighborhood size.
    n_pcs:
        Number of PCs to use. None uses all of X_pca.
    random_state:
        Seed for reproducibility.

    Returns
    -------
    AnnData with .uns['neighbors'] and .obsp graph matrices added in-place.
    """
    if "X_pca" not in adata.obsm:
        raise ValueError(
            "neighbors requires a PCA embedding (adata.obsm['X_pca']). "
            "Run the 'pca' primitive before 'neighbors'."
        )

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=random_state,
    )
    logger.info(f"neighbors: n_neighbors={n_neighbors}, n_pcs={n_pcs}, random_state={random_state}")
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("neighbors")
def _check_neighbors(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate neighbor-graph construction.

    Checks:
    1. n_neighbors is comparable to or exceeds the number of cells.
    2. connectivities matrix is missing (graph not built).
    """
    warnings: list[Warning_] = []

    n_neighbors = int(params.get("n_neighbors", 15))
    n_obs = adata_after.n_obs
    if n_neighbors >= n_obs:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"n_neighbors={n_neighbors} is >= the number of cells ({n_obs}). "
                    "The neighbor graph will be degenerate; reduce n_neighbors."
                ),
                primitive="neighbors",
            )
        )

    if "connectivities" not in adata_after.obsp:
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "Neighbor graph connectivities were not written to .obsp. "
                    "Downstream UMAP and clustering will fail."
                ),
                primitive="neighbors",
            )
        )

    return warnings

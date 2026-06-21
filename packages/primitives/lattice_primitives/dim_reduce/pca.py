"""Primitive: pca.

Wraps sc.pp.pca to compute a principal-component embedding restricted to the
highly variable genes. PCA is the standard linear dimensionality-reduction step
that feeds neighbor-graph construction and clustering.

Requires highly_variable_genes to have run first (uses mask_var="highly_variable").

Additions to AnnData:
    adata.obsm["X_pca"]    — n_obs x n_comps PC coordinates
    adata.varm["PCs"]      — gene loadings (documented; not schema-validated)
    adata.uns["pca"]       — {'variance', 'variance_ratio'}

Layer state is passthrough (log_normalized in, log_normalized out).
random_state is threaded for reproducibility (AGENTS.md: non-negotiable).
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
    var_columns=["highly_variable"],  # produced by highly_variable_genes
)
_PRODUCES = AnnDataState(
    adds_obsm=["X_pca"],
    adds_uns=["pca"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "n_comps": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 50,
            "description": (
                "Number of principal components. 50 is the scanpy/Theis convention. "
                "Automatically capped at min(n_obs, n_vars) - 1 for small datasets."
            ),
            "min_value": 2,
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
    name="pca",
    category="dim_reduce",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_pca",
    row_modifying=False,
)
def pca(
    adata: AnnData,
    n_comps: int = 50,
    random_state: int = 0,
) -> AnnData:
    """Compute PCA on highly variable genes via sc.pp.pca.

    Parameters
    ----------
    adata:
        AnnData in log_normalized layer_state with .var['highly_variable'] set.
    n_comps:
        Number of components. Capped at min(n_obs, n_vars) - 1 to stay valid on
        small datasets.
    random_state:
        Seed for the (randomized) SVD solver. Threaded for reproducibility.

    Returns
    -------
    AnnData with .obsm['X_pca'], .varm['PCs'], .uns['pca'] added in-place.
    """
    # Cap n_comps so PCA stays valid on small inputs (e.g. test fixtures).
    max_comps = min(adata.n_obs, adata.n_vars) - 1
    n_comps_effective = min(n_comps, max_comps)
    if n_comps_effective < n_comps:
        logger.info(
            f"pca: capped n_comps from {n_comps} to {n_comps_effective} "
            f"(dataset is {adata.n_obs}x{adata.n_vars})"
        )

    sc.pp.pca(
        adata,
        n_comps=n_comps_effective,
        mask_var="highly_variable",
        random_state=random_state,
    )
    logger.info(
        f"pca: computed {n_comps_effective} components on highly variable genes "
        f"(random_state={random_state})"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("pca")
def _check_pca(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate PCA output.

    Checks:
    1. NaNs in X_pca (degenerate input).
    2. PC1 explains a dominant fraction of variance (>0.5) — possible technical
       artifact / single overwhelming axis.
    """
    warnings: list[Warning_] = []

    x_pca = adata_after.obsm.get("X_pca")
    if x_pca is not None and not np.all(np.isfinite(np.asarray(x_pca))):
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "X_pca contains non-finite values (NaN/Inf). The input matrix "
                    "may contain NaNs or all-zero genes. Verify normalization."
                ),
                primitive="pca",
            )
        )

    variance_ratio = adata_after.uns.get("pca", {}).get("variance_ratio")
    if variance_ratio is not None and len(variance_ratio) > 0:
        pc1 = float(variance_ratio[0])
        if pc1 > 0.5:
            warnings.append(
                Warning_(
                    severity="warn",
                    message=(
                        f"PC1 explains {pc1 * 100:.0f}% of variance (>50%). A single "
                        "dominant axis can indicate a strong technical effect, an "
                        "outlier population, or insufficient feature selection."
                    ),
                    primitive="pca",
                )
            )

    return warnings

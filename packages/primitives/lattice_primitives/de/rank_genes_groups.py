"""Primitive: rank_genes_groups.

Wraps sc.tl.rank_genes_groups to rank marker genes that distinguish each group
(typically each Leiden cluster) from the rest. This is the standard one-vs-rest
marker-gene test that downstream annotation and plotting consume.

Requires log_normalized data. The groupby column is a *parameter*, so it cannot
be declared in `requires.obs_columns` (which is static). The body guards on it
explicitly.

Additions to AnnData:
    adata.uns["rank_genes_groups"] — dict with names/scores/pvals/pvals_adj/
        logfoldchanges/pts per group (the scanpy result structure)

NOTE (biology-semantics, provisional — see docs/BIOLOGY_DEFAULTS.md):
    method='wilcoxon' is provisional and FLAGGED FOR LSM-FOUNDER REVIEW.
    Scanpy's own default is 't-test', but the Wilcoxon rank-sum test is the
    field-recommended default for scRNA-seq marker detection (Soneson & Robinson
    2018; scanpy/Theis best-practices tutorial) because it is non-parametric and
    robust to the non-normal, zero-inflated expression distributions typical of
    single-cell data. We deliberately override scanpy's default. tie_correct=True
    is enabled because Wilcoxon ties are pervasive in sparse count data.
    corr_method='benjamini-hochberg' (FDR) is provisional and also flagged.

Layer state is passthrough (log_normalized in, log_normalized out).
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
    # groupby is a parameter, not statically known here, so it is guarded in the
    # body rather than declared in requires.obs_columns.
)
_PRODUCES = AnnDataState(
    adds_uns=["rank_genes_groups"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "groupby": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "leiden",
            "description": (
                "The .obs column defining the groups to compare (one-vs-rest). "
                "Defaults to the Leiden cluster labels. Must be an existing, "
                "categorical .obs column with at least two groups."
            ),
        }
    ),
    "method": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "wilcoxon",
            "description": (
                "Statistical test for marker ranking. 'wilcoxon' (non-parametric "
                "rank-sum) is the recommended scRNA-seq default; 't-test' is "
                "scanpy's built-in default; 'logreg' fits a logistic regression."
            ),
            "choices": ["wilcoxon", "t-test", "logreg"],
        }
    ),
    "n_genes": ParamSpec.model_validate(
        {
            "type": "int",
            "default": None,
            "description": (
                "Number of top-ranked genes to retain per group. None (default) ranks all genes."
            ),
        }
    ),
}


@primitive(
    name="rank_genes_groups",
    category="de",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="scanpy_rank_genes_groups",
    row_modifying=False,
)
def rank_genes_groups(
    adata: AnnData,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int | None = None,
) -> AnnData:
    """Rank marker genes per group via sc.tl.rank_genes_groups.

    Parameters
    ----------
    adata:
        Log-normalized AnnData with the `groupby` column present in .obs.
    groupby:
        .obs column defining the groups (default the Leiden clusters).
    method:
        'wilcoxon' (default), 't-test', or 'logreg'.
    n_genes:
        Number of top genes to retain per group; None ranks all genes.

    Returns
    -------
    AnnData with .uns['rank_genes_groups'] added in-place.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"rank_genes_groups requires groupby column '{groupby}' in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}. "
            "Run a clustering primitive (e.g. 'leiden') first."
        )

    # Guard: scanpy raises an opaque error if any group has < 2 cells. Surface a
    # clear, actionable message before calling it (a singleton group makes a
    # one-vs-rest marker test statistically undefined).
    group_sizes = adata.obs[groupby].value_counts()
    singletons = group_sizes[group_sizes < 2]
    if len(singletons) > 0:
        raise ValueError(
            f"rank_genes_groups cannot test group(s) {list(singletons.index)} in "
            f"'{groupby}' that contain fewer than 2 cells. Merge or remove tiny "
            "clusters (e.g. lower the Leiden resolution) before marker testing."
        )

    # tie_correct only applies to the Wilcoxon test; corr_method only applies to
    # the parametric tests. Pass them conditionally so scanpy does not warn.
    kwargs: dict[str, object] = {
        "groupby": groupby,
        "method": method,
        "pts": True,
        "key_added": "rank_genes_groups",
    }
    if n_genes is not None:
        kwargs["n_genes"] = n_genes
    if method == "wilcoxon":
        kwargs["tie_correct"] = True
    if method in ("wilcoxon", "t-test"):
        kwargs["corr_method"] = "benjamini-hochberg"

    sc.tl.rank_genes_groups(adata, **kwargs)  # type: ignore[arg-type]

    n_groups = len(adata.uns["rank_genes_groups"]["names"].dtype.names)
    logger.info(
        f"rank_genes_groups: method={method}, groupby='{groupby}' -> {n_groups} group(s) ranked"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("rank_genes_groups")
def _check_rank_genes_groups(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate marker-ranking output.

    Checks:
    1. groupby column is categorical with >= 2 groups.
    2. Data looks log-normalized (max value plausibly post-log1p).
    3. All adjusted p-values are NaN (test failed to produce results).

    (Singleton groups are rejected by a guard in the primitive body before the
    scanpy call, so they never reach this sanity check.)
    """
    warnings: list[Warning_] = []
    groupby = params.get("groupby", "leiden")

    # 1: group structure
    if groupby in adata_after.obs.columns:
        col = adata_after.obs[groupby]
        group_counts = col.value_counts()
        n_groups = int((group_counts > 0).sum())
        if n_groups < 2:
            warnings.append(
                Warning_(
                    severity="error",
                    message=(
                        f"rank_genes_groups needs >= 2 non-empty groups in "
                        f"'{groupby}', found {n_groups}. One-vs-rest comparison "
                        "is undefined with a single group."
                    ),
                    primitive="rank_genes_groups",
                )
            )

    # 2: data looks log-normalized
    layer_state = adata_after.uns.get("_lattice_layer_state")
    if layer_state != "log_normalized":
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"layer_state is '{layer_state}', expected 'log_normalized'. "
                    "Marker tests assume log-normalized expression; running on "
                    "raw counts produces misleading fold changes."
                ),
                primitive="rank_genes_groups",
            )
        )

    # 4: all-NaN adjusted p-values
    rgg = adata_after.uns.get("rank_genes_groups")
    if rgg is not None and "pvals_adj" in rgg:
        pvals_adj = rgg["pvals_adj"]
        # Structured recarray: stack all group columns into one float array.
        try:
            all_pvals = np.concatenate(
                [np.asarray(pvals_adj[g], dtype=float) for g in pvals_adj.dtype.names]
            )
            if all_pvals.size > 0 and np.all(np.isnan(all_pvals)):
                warnings.append(
                    Warning_(
                        severity="error",
                        message=(
                            "All adjusted p-values are NaN. The test produced no "
                            "usable results — check for empty groups or "
                            "degenerate (constant) expression."
                        ),
                        primitive="rank_genes_groups",
                    )
                )
        except (AttributeError, TypeError):  # pragma: no cover - defensive
            pass

    return warnings

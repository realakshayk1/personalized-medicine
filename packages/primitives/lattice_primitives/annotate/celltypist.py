"""Primitive: celltypist.

Wraps `celltypist.annotate(adata, model=..., majority_voting=True)` to assign a
cell-type label to every cell using a pre-trained CellTypist logistic-regression
model, then writes the predicted labels back onto adata.obs.

NOTE (biology-semantics, PROVISIONAL — see docs/BIOLOGY_DEFAULTS.md §8):
    The default model "Immune_All_Low.pkl" is PROVISIONAL and FLAGGED FOR
    LSM-FOUNDER REVIEW. It is an immune-focused model; tissue-specific models may
    be more appropriate for non-immune data. majority_voting=True (PROVISIONAL)
    refines per-cell predictions over over-clustered neighborhoods.

HARD INPUT REQUIREMENT (PROVISIONAL — flagged for LSM review):
    CellTypist expects expression log1p-normalized to a fixed target_sum of 1e4
    ("log1p-CP10K"). This DIFFERS from Lattice's default median normalization
    (target_sum=None). Re-normalizing adata.X in place would corrupt the lineage
    tag and break every downstream primitive that assumes the median-normalized
    state. To stay safe we build a 1e4-normalized COPY locally for prediction
    ONLY, and never mutate adata.X. We require the incoming layer_state to be
    "log_normalized" (we re-derive counts via expm1 from the log data, then
    re-normalize the copy to 1e4). This re-derivation is approximate; a dedicated
    raw-counts layer would be cleaner. FLAGGED for LSM review.

Columns added to .obs:
    celltypist_labels            — categorical predicted cell-type label per cell
                                   (majority-voted when majority_voting=True)
    celltypist_conf_score        — float prediction confidence/probability per cell

adata.uns additions:
    celltypist                   — dict with model name, majority_voting flag, and
                                   the set of predicted labels (summary only, no
                                   raw matrices — avoids the context-bomb failure)

Layer state is passthrough (the original adata.X is untouched).
"""

from __future__ import annotations

from typing import Any

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
    layer_state="log_normalized",
)
_PRODUCES = AnnDataState(
    obs_columns=["celltypist_labels", "celltypist_conf_score"],
    adds_uns=["celltypist"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "model": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "Immune_All_Low.pkl",
            "description": (
                "CellTypist model name (.pkl) or path. PROVISIONAL default "
                "'Immune_All_Low.pkl' (immune-focused) — flagged for LSM review; "
                "tissue-specific models may suit non-immune data better."
            ),
        }
    ),
    "majority_voting": ParamSpec.model_validate(
        {
            "type": "bool",
            "default": True,
            "description": (
                "If True, refine per-cell predictions by majority vote over "
                "over-clustered neighborhoods (CellTypist majority_voting). "
                "PROVISIONAL default — flagged for LSM review."
            ),
        }
    ),
}


def _to_log1p_cp10k(adata: AnnData) -> AnnData:
    """Return a COPY of adata re-normalized to log1p-CP10K (target_sum=1e4).

    The incoming adata is log_normalized (median target). We invert log1p with
    expm1 to recover relative counts, then renormalize to 1e4 and re-log. We
    never touch the input adata's .X. Only X is needed for prediction.
    """
    work = AnnData(
        X=adata.X.copy(),
        obs=adata.obs[[]].copy(),  # keep index only; CellTypist needs obs_names
        var=adata.var[[]].copy(),  # keep index only; CellTypist matches var_names
    )
    work.obs_names = adata.obs_names
    work.var_names = adata.var_names

    # Invert the log1p that produced log_normalized data, then renormalize.
    if scipy.sparse.issparse(work.X):
        work.X = work.X.copy()
        work.X.data = np.expm1(work.X.data)
    else:
        work.X = np.expm1(np.asarray(work.X))
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    return work


@primitive(
    name="celltypist",
    category="annotate",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="dominguez_celltypist_2022",
    row_modifying=False,
)
def celltypist(
    adata: AnnData,
    model: str | Any = "Immune_All_Low.pkl",
    majority_voting: bool = True,
) -> AnnData:
    """Annotate cell types with CellTypist and write labels back to adata.obs.

    Parameters
    ----------
    adata:
        Log-normalized AnnData (layer_state='log_normalized'). A 1e4-normalized
        COPY is built internally for prediction; adata.X is NOT modified.
    model:
        CellTypist model name (.pkl) or filesystem path to a model. PROVISIONAL
        default 'Immune_All_Low.pkl'.
    majority_voting:
        Refine predictions over neighborhoods (PROVISIONAL default True).

    Returns
    -------
    AnnData with .obs['celltypist_labels'], .obs['celltypist_conf_score'] and
    .uns['celltypist'] added in-place. The original .X is untouched.

    Notes
    -----
    Importing celltypist and loading a model may require a network download on
    first use. The `model` argument also accepts a local filesystem path so the
    primitive can run fully offline once a model is available locally.
    """
    import celltypist
    from celltypist import models as ct_models

    # Build the 1e4-normalized copy CellTypist expects (does not mutate adata.X).
    work = _to_log1p_cp10k(adata)

    # Resolve the model: accept either an already-loaded Model (useful for tests
    # and warm sessions), or a name/path that celltypist loads (downloading if
    # needed). This keeps the primitive runnable offline given a local model.
    if isinstance(model, ct_models.Model):
        loaded_model = model
    else:
        loaded_model = ct_models.Model.load(model=model)

    predictions = celltypist.annotate(
        work,
        model=loaded_model,
        majority_voting=majority_voting,
    )

    pred_obs = predictions.predicted_labels
    # When majority_voting=True, predicted_labels is a DataFrame with several
    # columns ('predicted_labels', 'over_clustering', 'majority_voting'); use the
    # majority-voted column. Otherwise it is a single-column DataFrame/Series.
    if hasattr(pred_obs, "columns") and "majority_voting" in getattr(pred_obs, "columns", []):
        labels = pred_obs["majority_voting"]
    elif hasattr(pred_obs, "columns") and "predicted_labels" in getattr(pred_obs, "columns", []):
        labels = pred_obs["predicted_labels"]
    else:
        labels = pred_obs

    labels = np.asarray(labels).reshape(-1).astype(str)
    adata.obs["celltypist_labels"] = labels
    adata.obs["celltypist_labels"] = adata.obs["celltypist_labels"].astype("category")

    # Per-cell confidence: max class probability from the decision matrix.
    prob = predictions.probability_matrix
    try:
        conf = np.asarray(prob).max(axis=1).reshape(-1).astype(float)
    except Exception:  # pragma: no cover - defensive
        conf = np.full(adata.n_obs, np.nan, dtype=float)
    adata.obs["celltypist_conf_score"] = conf

    unique_labels = sorted(set(labels.tolist()))
    adata.uns["celltypist"] = {
        "model": str(model),
        "majority_voting": bool(majority_voting),
        "labels": unique_labels,
        "n_labels": len(unique_labels),
    }
    logger.info(
        f"celltypist: model={model}, majority_voting={majority_voting} -> "
        f"{len(unique_labels)} distinct label(s) across {adata.n_obs} cells"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("celltypist")
def _check_celltypist(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate CellTypist annotation output.

    Checks:
    1. Labels column present and non-empty (error if missing).
    2. A single label dominates (>95% of cells) — likely wrong model for the
       tissue (warn).
    3. Low mean confidence (< 0.5) — predictions are weak (warn).
    """
    warnings: list[Warning_] = []

    if "celltypist_labels" not in adata_after.obs.columns:
        warnings.append(
            Warning_(
                severity="error",
                message="CellTypist produced no 'celltypist_labels' column.",
                primitive="celltypist",
            )
        )
        return warnings

    labels = adata_after.obs["celltypist_labels"]
    if labels.isna().all():
        warnings.append(
            Warning_(
                severity="error",
                message="All CellTypist labels are missing (annotation failed).",
                primitive="celltypist",
            )
        )
        return warnings

    largest_frac = float(labels.value_counts(normalize=True).iloc[0])
    if largest_frac > 0.95:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"A single cell type holds {largest_frac * 100:.0f}% of cells. "
                    "The chosen CellTypist model may not match this tissue; "
                    "consider a tissue-specific model (flagged for LSM review)."
                ),
                primitive="celltypist",
            )
        )

    if "celltypist_conf_score" in adata_after.obs.columns:
        conf = adata_after.obs["celltypist_conf_score"].to_numpy(dtype=float)
        if conf.size > 0 and np.isfinite(conf).any():
            mean_conf = float(np.nanmean(conf))
            if mean_conf < 0.5:
                warnings.append(
                    Warning_(
                        severity="warn",
                        message=(
                            f"Mean CellTypist confidence is {mean_conf:.2f} (< 0.50). "
                            "Predictions are weak; review labels before trusting them."
                        ),
                        primitive="celltypist",
                    )
                )

    return warnings

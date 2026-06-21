"""Primitive: mllm_consensus  (the mLLMCelltype consensus pattern).

Annotates each cluster's cell type by asking TWO language models independently
(Opus + Sonnet) to name the cell type from its top marker genes, then computing a
consensus across their answers. This is the ONE sanctioned LLM-driven primitive
(AGENTS.md): the models pick a NAME from cited markers — they never write or run
analysis code (ADR-001). Output still flows through schema validation and a
sanity check like every other primitive.

Inputs:
    - adata.uns['rank_genes_groups'] (run the 'rank_genes_groups' primitive
      first). The top-N marker genes per group are extracted and sent to the
      models. The body guards on this explicitly.
    - groupby: the .obs column whose clusters were ranked (must exist).

Algorithm (mLLMCelltype, Yang et al. 2025):
    For each cluster:
      1. Build a structured prompt (top markers + species + tissue context).
      2. Query two annotators independently -> two predictions.
      3. Consensus proportion (CP) = fraction of annotators agreeing on the
         plurality label.
      4. Shannon entropy H over the prediction distribution (uncertainty).
      5. Confidence bucket:  CP > 0.8 -> "high"; 0.5 <= CP <= 0.8 -> "review";
         CP < 0.5 -> "needs-expert".
      6. Markers cited in the prediction are recorded (required).

OPEN GAP (PROVISIONAL — flagged prominently for LSM review, see
docs/BIOLOGY_DEFAULTS.md §8): mLLMCelltype emits FREE-TEXT cell-type names, NOT
Cell Ontology IDs (CL:XXXXXXX). AGENTS.md failure-mode #4 requires a CL ID (with
cited markers) or the label MUST be flagged "needs review". We provide a
name->CL mapping HOOK (`cl_map`, an empty dict by default = no-op). When no CL ID
resolves for a cluster, its confidence is downgraded to "needs-expert"
regardless of CP, and a flag is recorded. A real CellOntologyMapper is a
deliberate follow-up.

CRITICAL DESIGN (testability, no API key): the LLM call is INJECTABLE. The
`annotator` parameter is a callable(cluster_context: dict) -> dict that defaults
to a real Claude-backed implementation via claude-agent-sdk if importable, else a
stub that raises a clear error. Tests inject a deterministic fake annotator so the
consensus/entropy/confidence logic is fully covered with no network and no
ANTHROPIC_API_KEY.

Columns added to .obs:
    cell_type              — categorical consensus label per cell (cluster -> label)
    cell_type_confidence   — categorical: "high" / "review" / "needs-expert" per cell

adata.uns additions:
    mllm_consensus         — dict: per-cluster reasoning, markers, predictions,
                             CP, entropy, confidence, resolved CL id (summary only;
                             no raw matrices — avoids the context-bomb failure)

Layer state is passthrough; adata.X is untouched.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable
from typing import Any

import numpy as np
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_
from lattice_primitives.registry import primitive, sanity_check

# ---------------------------------------------------------------------------
# Annotator: injectable LLM call
# ---------------------------------------------------------------------------

# A prediction is {"label": str, "cited_markers": list[str], "model": str,
#                  "reasoning": str}. An annotator takes a cluster context dict
# {"cluster": str, "markers": list[str], "species": str, "tissue": str} and
# returns ONE prediction. The default annotator queries two Claude models.

PredictionDict = dict[str, Any]
AnnotatorFn = Callable[[dict[str, Any]], list[PredictionDict]]


def _default_claude_annotator(context: dict[str, Any]) -> list[PredictionDict]:
    """Query Opus and Sonnet independently for a cluster's cell type.

    Defaults to a real claude-agent-sdk-backed call when the SDK is importable.
    If it is not importable (or no API key is configured) this raises a clear
    error — callers running offline MUST inject their own annotator.

    NOTE: this is the sanctioned LLM step. The model returns only a NAME and the
    markers it relied on; it never emits analysis code.
    """
    try:
        import claude_agent_sdk  # noqa: F401
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "No annotator configured: claude-agent-sdk is not importable and no "
            "`annotator` callable was injected. Provide annotator=<fn> to run "
            "mllm_consensus offline, or install claude-agent-sdk."
        ) from exc

    # Real implementation would issue two independent SDK queries (Opus, Sonnet)
    # with the structured prompt built from `context`, parse each into a
    # PredictionDict, and return both. Kept as an explicit gap so the primitive
    # never silently fabricates predictions without a configured backend.
    raise RuntimeError(
        "No annotator configured: the default Claude-backed annotator requires a "
        "live claude-agent-sdk session and ANTHROPIC_API_KEY. Inject an "
        "`annotator` callable for tests or offline runs."
    )


def _build_cluster_context(
    cluster: str,
    markers: list[str],
    species: str,
    tissue: str,
) -> dict[str, Any]:
    """Assemble the structured per-cluster prompt context."""
    return {
        "cluster": str(cluster),
        "markers": list(markers),
        "species": species,
        "tissue": tissue,
        "prompt": (
            f"Species: {species}. Tissue: {tissue}. "
            f"Cluster {cluster} top marker genes (ranked): {', '.join(markers)}. "
            "Name the single most likely cell type. Cite the specific marker "
            "genes that justify the call. Do not write code."
        ),
    }


def _consensus(predictions: list[PredictionDict]) -> tuple[str, float, float, list[str]]:
    """Compute (label, consensus_proportion, shannon_entropy, cited_markers).

    CP = count(plurality label) / n_predictions.
    Entropy is computed over the label distribution (natural log, in nats).
    Cited markers are the union of markers cited by annotators that voted for the
    winning label.
    """
    labels = [str(p.get("label", "")).strip() for p in predictions]
    n = len(labels)
    counts = Counter(labels)
    winner, winner_count = counts.most_common(1)[0]
    cp = winner_count / n if n else 0.0

    # Shannon entropy over the label distribution.
    entropy = 0.0
    for _, c in counts.items():
        p = c / n
        entropy -= p * math.log(p)

    cited: list[str] = []
    for pred in predictions:
        if str(pred.get("label", "")).strip() == winner:
            for m in pred.get("cited_markers", []) or []:
                if m not in cited:
                    cited.append(str(m))
    return winner, cp, entropy, cited


def _confidence_bucket(cp: float) -> str:
    """Map a consensus proportion to a confidence bucket."""
    if cp > 0.8:
        return "high"
    if cp >= 0.5:
        return "review"
    return "needs-expert"


# ---------------------------------------------------------------------------
# Primitive registration
# ---------------------------------------------------------------------------

_REQUIRES = AnnDataState(
    layer_state="log_normalized",
    # rank_genes_groups + groupby are guarded in the body (uns / parameter).
)
_PRODUCES = AnnDataState(
    obs_columns=["cell_type", "cell_type_confidence"],
    adds_uns=["mllm_consensus"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "groupby": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "leiden",
            "description": (
                "The .obs cluster column that was marker-ranked. Must match the "
                "groupby used for rank_genes_groups."
            ),
        }
    ),
    "n_top_markers": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 25,
            "description": "Number of top marker genes per cluster sent to the models.",
            "min_value": 1,
        }
    ),
    "species": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "human",
            "description": "Species context for the annotation prompt. PROVISIONAL.",
        }
    ),
    "tissue": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "unknown",
            "description": "Tissue context for the annotation prompt. PROVISIONAL.",
        }
    ),
}


@primitive(
    name="mllm_consensus",
    category="annotate",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="yang_mllmcelltype_2025",
    row_modifying=False,
)
def mllm_consensus(
    adata: AnnData,
    groupby: str = "leiden",
    n_top_markers: int = 25,
    species: str = "human",
    tissue: str = "unknown",
    annotator: AnnotatorFn | None = None,
    cl_map: dict[str, str] | None = None,
) -> AnnData:
    """Annotate clusters via multi-LLM consensus over cited markers.

    Parameters
    ----------
    adata:
        Log-normalized AnnData with .uns['rank_genes_groups'] present.
    groupby:
        Cluster column that was marker-ranked.
    n_top_markers:
        Top markers per cluster sent to the annotators.
    species, tissue:
        Context strings injected into the prompt (PROVISIONAL defaults).
    annotator:
        Injectable callable(context) -> list[prediction]. Defaults to a
        Claude-backed annotator (requires claude-agent-sdk + API key). Tests
        inject a deterministic fake.
    cl_map:
        Name -> Cell Ontology ID mapping hook. Empty/None = no-op (the OPEN GAP).
        Clusters with no resolved CL id are forced to "needs-expert".

    Returns
    -------
    AnnData with .obs['cell_type'], .obs['cell_type_confidence'] and
    .uns['mllm_consensus'] added in-place.
    """
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "mllm_consensus requires marker genes in adata.uns['rank_genes_groups']. "
            "Run the 'rank_genes_groups' primitive before 'mllm_consensus'."
        )
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"mllm_consensus requires groupby column '{groupby}' in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}."
        )

    if annotator is None:
        annotator = _default_claude_annotator
    cl_map = cl_map or {}

    rgg = adata.uns["rank_genes_groups"]
    names = rgg["names"]  # structured recarray: one column per group
    group_keys = list(names.dtype.names)

    cluster_records: dict[str, dict[str, Any]] = {}
    label_per_cluster: dict[str, str] = {}
    conf_per_cluster: dict[str, str] = {}

    for grp in group_keys:
        top_markers = [str(g) for g in np.asarray(names[grp])[:n_top_markers]]
        context = _build_cluster_context(grp, top_markers, species, tissue)
        predictions = annotator(context)
        if not predictions:
            raise ValueError(f"Annotator returned no predictions for cluster '{grp}'.")

        label, cp, entropy, cited = _consensus(predictions)
        confidence = _confidence_bucket(cp)

        # OPEN GAP: resolve a Cell Ontology id, else force needs-expert.
        cl_id = cl_map.get(label)
        cl_resolved = cl_id is not None
        if not cl_resolved:
            confidence = "needs-expert"

        label_per_cluster[grp] = label
        conf_per_cluster[grp] = confidence
        cluster_records[grp] = {
            "label": label,
            "consensus_proportion": cp,
            "shannon_entropy": entropy,
            "confidence": confidence,
            "cited_markers": cited,
            "top_markers": top_markers,
            "cl_id": cl_id,
            "cl_resolved": cl_resolved,
            "predictions": [
                {
                    "model": str(p.get("model", "unknown")),
                    "label": str(p.get("label", "")),
                    "cited_markers": list(p.get("cited_markers", []) or []),
                    "reasoning": str(p.get("reasoning", "")),
                }
                for p in predictions
            ],
        }

    # Map cluster-level results onto each cell.
    grp_values = adata.obs[groupby].astype(str)
    adata.obs["cell_type"] = grp_values.map(label_per_cluster).astype("category")
    adata.obs["cell_type_confidence"] = grp_values.map(conf_per_cluster).astype("category")

    adata.uns["mllm_consensus"] = {
        "groupby": groupby,
        "species": species,
        "tissue": tissue,
        "n_top_markers": int(n_top_markers),
        "clusters": cluster_records,
        "cl_mapping_used": bool(cl_map),
    }

    n_high = sum(1 for v in conf_per_cluster.values() if v == "high")
    logger.info(
        f"mllm_consensus: groupby='{groupby}' -> {len(group_keys)} cluster(s) "
        f"annotated ({n_high} high-confidence). CL mapping "
        f"{'present' if cl_map else 'EMPTY (open gap; all forced needs-expert)'}."
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("mllm_consensus")
def _check_mllm_consensus(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate consensus annotation output.

    Checks:
    1. No clusters annotated at all (error).
    2. Any cluster with confidence "needs-expert" or "review", OR missing a
       Cell Ontology id (warn — must be surfaced for review per failure-mode #4).
    """
    warnings: list[Warning_] = []

    record = adata_after.uns.get("mllm_consensus")
    clusters = (record or {}).get("clusters", {})
    if not clusters:
        warnings.append(
            Warning_(
                severity="error",
                message="mllm_consensus annotated no clusters (empty result).",
                primitive="mllm_consensus",
            )
        )
        return warnings

    flagged: list[str] = []
    missing_cl: list[str] = []
    for grp, rec in clusters.items():
        if rec.get("confidence") in ("needs-expert", "review"):
            flagged.append(str(grp))
        if not rec.get("cl_resolved", False):
            missing_cl.append(str(grp))

    if flagged:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Cluster(s) {sorted(set(flagged))} have low-confidence "
                    "consensus annotations (review / needs-expert). Confirm these "
                    "cell types before reporting."
                ),
                primitive="mllm_consensus",
            )
        )
    if missing_cl:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Cluster(s) {sorted(set(missing_cl))} have no resolved Cell "
                    "Ontology ID (CL:XXXXXXX). Per the annotation contract these "
                    "labels are 'needs review' until a CL id is assigned "
                    "(open gap: name->CL mapping)."
                ),
                primitive="mllm_consensus",
            )
        )

    return warnings

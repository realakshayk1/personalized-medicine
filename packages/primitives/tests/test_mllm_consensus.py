"""Tests for the mllm_consensus (multi-LLM consensus annotation) primitive.

A deterministic FAKE annotator is injected so the consensus / entropy /
confidence logic is fully exercised with no network and no ANTHROPIC_API_KEY.
"""

from __future__ import annotations

from typing import Any

import pytest
from anndata import AnnData

import lattice_primitives.annotate.mllm_consensus  # noqa: F401  (register)
import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.de.rank_genes_groups  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
from lattice_primitives import run_primitive
from lattice_primitives.annotate.mllm_consensus import (
    _confidence_bucket,
    _consensus,
)

# ---------------------------------------------------------------------------
# Pipeline helper: cluster + rank markers (local to this test file)
# ---------------------------------------------------------------------------


def _through_ranked(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    out, _, _ = run_primitive("rank_genes_groups", out, groupby="leiden", method="wilcoxon")
    return out


def _agree_annotator(context: dict[str, Any]) -> list[dict[str, Any]]:
    """Two models that always agree -> CP=1.0, entropy=0, confidence high (if CL)."""
    label = f"cell_type_{context['cluster']}"
    markers = context["markers"][:3]
    return [
        {"model": "opus", "label": label, "cited_markers": markers, "reasoning": "r1"},
        {"model": "sonnet", "label": label, "cited_markers": markers, "reasoning": "r2"},
    ]


def _disagree_annotator(context: dict[str, Any]) -> list[dict[str, Any]]:
    """Two models that always disagree -> CP=0.5, entropy>0, confidence 'review'."""
    markers = context["markers"][:3]
    return [
        {"model": "opus", "label": "T cell", "cited_markers": markers, "reasoning": "a"},
        {"model": "sonnet", "label": "B cell", "cited_markers": markers, "reasoning": "b"},
    ]


# ---------------------------------------------------------------------------
# Consensus math (pure-function unit tests)
# ---------------------------------------------------------------------------


def test_consensus_full_agreement() -> None:
    preds = [
        {"label": "T cell", "cited_markers": ["CD3D", "CD3E"]},
        {"label": "T cell", "cited_markers": ["CD3D"]},
    ]
    label, cp, entropy, cited = _consensus(preds)
    assert label == "T cell"
    assert cp == 1.0
    assert entropy == pytest.approx(0.0)
    assert "CD3D" in cited and "CD3E" in cited


def test_consensus_split_vote() -> None:
    preds = [
        {"label": "T cell", "cited_markers": ["CD3D"]},
        {"label": "B cell", "cited_markers": ["MS4A1"]},
    ]
    label, cp, entropy, cited = _consensus(preds)
    assert cp == 0.5
    assert entropy == pytest.approx(0.6931, abs=1e-3)  # ln(2)
    # Cited markers come only from the winning label's voters.
    winner_markers = {"CD3D"} if label == "T cell" else {"MS4A1"}
    assert set(cited) == winner_markers


def test_confidence_buckets() -> None:
    assert _confidence_bucket(0.9) == "high"
    assert _confidence_bucket(0.81) == "high"
    assert _confidence_bucket(0.8) == "review"
    assert _confidence_bucket(0.5) == "review"
    assert _confidence_bucket(0.49) == "needs-expert"


# ---------------------------------------------------------------------------
# Primitive: happy path, guards, immutability, sanity triggers
# ---------------------------------------------------------------------------


def test_mllm_happy_path_with_cl_map(adata_log_normalized: AnnData) -> None:
    adata = _through_ranked(adata_log_normalized)
    groups = list(adata.uns["rank_genes_groups"]["names"].dtype.names)
    # Provide a CL mapping for every agreeing label so confidence can reach "high".
    cl_map = {f"cell_type_{g}": f"CL:000{i:04d}" for i, g in enumerate(groups)}

    adata_out, warnings, result = run_primitive(
        "mllm_consensus",
        adata,
        groupby="leiden",
        annotator=_agree_annotator,
        cl_map=cl_map,
    )
    assert "cell_type" in adata_out.obs.columns
    assert "cell_type_confidence" in adata_out.obs.columns
    assert "mllm_consensus" in adata_out.uns
    clusters = adata_out.uns["mllm_consensus"]["clusters"]
    # Full agreement + resolved CL id -> all high confidence, no warnings.
    assert all(rec["confidence"] == "high" for rec in clusters.values())
    assert all(rec["consensus_proportion"] == 1.0 for rec in clusters.values())
    assert all(rec["cited_markers"] for rec in clusters.values())  # markers cited
    assert all(w.severity != "error" for w in warnings)
    assert result.primitive == "mllm_consensus"


def test_mllm_missing_cl_forces_needs_expert(adata_log_normalized: AnnData) -> None:
    """Agreement is perfect, but with no CL map every cluster is needs-expert."""
    adata = _through_ranked(adata_log_normalized)
    adata_out, warnings, _ = run_primitive(
        "mllm_consensus", adata, groupby="leiden", annotator=_agree_annotator
    )
    clusters = adata_out.uns["mllm_consensus"]["clusters"]
    assert all(rec["confidence"] == "needs-expert" for rec in clusters.values())
    assert all(not rec["cl_resolved"] for rec in clusters.values())
    # Sanity check must surface low-confidence + missing-CL warnings.
    msgs = " ".join(w.message for w in warnings)
    assert any(w.severity == "warn" for w in warnings)
    assert "Cell Ontology" in msgs


def test_mllm_disagreement_is_review(adata_log_normalized: AnnData) -> None:
    """Split vote -> CP 0.5 -> 'review' bucket (when a CL id resolves)."""
    adata = _through_ranked(adata_log_normalized)
    cl_map = {"T cell": "CL:0000084", "B cell": "CL:0000236"}
    adata_out, warnings, _ = run_primitive(
        "mllm_consensus",
        adata,
        groupby="leiden",
        annotator=_disagree_annotator,
        cl_map=cl_map,
    )
    clusters = adata_out.uns["mllm_consensus"]["clusters"]
    assert all(rec["consensus_proportion"] == 0.5 for rec in clusters.values())
    assert all(rec["confidence"] == "review" for rec in clusters.values())
    assert all(rec["shannon_entropy"] > 0 for rec in clusters.values())
    assert any(w.severity == "warn" for w in warnings)


def test_mllm_rejects_missing_rank_genes_groups(adata_log_normalized: AnnData) -> None:
    """No rank_genes_groups in uns -> clear ValueError guard."""
    # leiden present but markers never ranked.
    out, _, _ = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    with pytest.raises(ValueError, match="rank_genes_groups"):
        run_primitive("mllm_consensus", out, groupby="leiden", annotator=_agree_annotator)


def test_mllm_preserves_condition(adata_log_normalized: AnnData) -> None:
    adata = _through_ranked(adata_log_normalized)
    before = list(adata.obs["condition"])
    out, _, _ = run_primitive("mllm_consensus", adata, groupby="leiden", annotator=_agree_annotator)
    assert list(out.obs["condition"]) == before

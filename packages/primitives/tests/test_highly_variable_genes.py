"""Tests for the highly_variable_genes primitive."""

from __future__ import annotations

import pytest
from anndata import AnnData

import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
from lattice_primitives import run_primitive
from lattice_primitives.registry import SchemaViolation


def test_happy_path_flags_genes(adata_log_normalized: AnnData) -> None:
    adata_out, warnings, result = run_primitive(
        "highly_variable_genes", adata_log_normalized, n_top_genes=200
    )
    assert "highly_variable" in adata_out.var.columns
    assert "dispersions_norm" in adata_out.var.columns
    assert int(adata_out.var["highly_variable"].sum()) == 200
    # layer_state is passthrough
    assert adata_out.uns["_lattice_layer_state"] == "log_normalized"
    assert result.primitive == "highly_variable_genes"


def test_rejects_wrong_layer_state(adata_raw_counts_state: AnnData) -> None:
    # requires log_normalized; raw_counts must be rejected before execution
    with pytest.raises(SchemaViolation):
        run_primitive("highly_variable_genes", adata_raw_counts_state, n_top_genes=50)


def test_sanity_warns_when_n_top_exceeds_genes(adata_log_normalized: AnnData) -> None:
    n_genes = adata_log_normalized.n_vars
    _, warnings, _ = run_primitive(
        "highly_variable_genes", adata_log_normalized, n_top_genes=n_genes + 1000
    )
    assert any("only" in w.message.lower() or "all genes" in w.message.lower() for w in warnings)


def test_explanation_rendered(adata_log_normalized: AnnData) -> None:
    _, _, result = run_primitive("highly_variable_genes", adata_log_normalized, n_top_genes=200)
    assert result.user_explanation
    assert "satija_hvg_2015" in result.methods_paragraph

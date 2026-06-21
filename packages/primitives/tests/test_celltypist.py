"""Tests for the celltypist annotation primitive.

CellTypist model loading needs a network download on first use. If a model is
not available offline, these tests SKIP gracefully (they must never hard-fail the
suite on a network error).

The synthetic AnnData here uses REAL gene symbols pulled from the loaded model's
feature list so that CellTypist has overlapping features to score (CellTypist
errors out if zero genes overlap the model). Data is constructed locally rather
than in conftest.py to keep these fixtures scoped to this file.
"""

from __future__ import annotations

import numpy as np
import pytest
import scanpy as sc
import scipy.sparse
from anndata import AnnData

import lattice_primitives.annotate.celltypist  # noqa: F401  (register)
import lattice_primitives.cluster.leiden  # noqa: F401  (register)
import lattice_primitives.dim_reduce.highly_variable_genes  # noqa: F401  (register)
import lattice_primitives.dim_reduce.neighbors  # noqa: F401  (register)
import lattice_primitives.dim_reduce.pca  # noqa: F401  (register)
from lattice_primitives import run_primitive

_MODEL_NAME = "Immune_All_Low.pkl"


def _load_model_or_skip():  # type: ignore[no-untyped-def]
    """Return a loaded CellTypist Model, or pytest.skip if unavailable offline."""
    try:
        import celltypist
        from celltypist import models as ct_models
    except Exception as exc:  # pragma: no cover - import guard
        pytest.skip(f"celltypist not importable: {exc}")

    # Try to load; if missing, attempt one download. Any failure -> skip.
    try:
        return ct_models.Model.load(model=_MODEL_NAME)
    except Exception:
        pass
    try:
        celltypist.models.download_models(model=_MODEL_NAME, force_update=False)
        return ct_models.Model.load(model=_MODEL_NAME)
    except Exception as exc:
        pytest.skip(
            f"CellTypist model '{_MODEL_NAME}' unavailable offline ({exc}); "
            "skipping celltypist annotation tests."
        )


def _log_normalized_with_model_genes(model, n_obs: int = 200, seed: int = 0) -> AnnData:  # type: ignore[no-untyped-def]
    """Build a log-normalized AnnData whose var_names are real model gene symbols.

    Uses the first ~400 features the model knows so CellTypist has features to
    score. Two synthetic populations give clustering real structure.
    """
    genes = [str(g) for g in list(model.features)[:400]]
    n_genes = len(genes)
    rng = np.random.default_rng(seed)
    n_a = n_obs // 2
    counts = rng.poisson(lam=3, size=(n_obs, n_genes)).astype(np.float32)
    counts[:n_a, : n_genes // 8] += rng.poisson(lam=20, size=(n_a, n_genes // 8)).astype(np.float32)
    counts[n_a:, n_genes // 8 : n_genes // 4] += rng.poisson(
        lam=20, size=(n_obs - n_a, n_genes // 8)
    ).astype(np.float32)

    adata = AnnData(X=scipy.sparse.csr_matrix(counts))
    adata.obs_names = [f"CELL_{i:04d}" for i in range(n_obs)]
    adata.var_names = genes
    adata.obs["condition"] = ["treated"] * n_a + ["control"] * (n_obs - n_a)
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)
    adata.uns["_lattice_layer_state"] = "log_normalized"
    return adata


def _clustered(adata: AnnData) -> AnnData:
    out, _, _ = run_primitive("highly_variable_genes", adata, n_top_genes=200)
    out, _, _ = run_primitive("pca", out, n_comps=30, random_state=0)
    out, _, _ = run_primitive("neighbors", out, n_neighbors=15, random_state=0)
    out, _, _ = run_primitive("leiden", out, resolution=1.0, random_state=0)
    return out


def test_celltypist_happy_path() -> None:
    model = _load_model_or_skip()
    adata = _clustered(_log_normalized_with_model_genes(model))
    x_before = adata.X.copy()

    adata_out, warnings, result = run_primitive(
        "celltypist", adata, model=model, majority_voting=True
    )
    assert "celltypist_labels" in adata_out.obs.columns
    assert "celltypist_conf_score" in adata_out.obs.columns
    assert "celltypist" in adata_out.uns
    assert adata_out.obs["celltypist_labels"].notna().any()
    assert result.primitive == "celltypist"
    assert all(w.severity != "error" for w in warnings)

    # adata.X must NOT be mutated (1e4 normalization happens on a local copy).
    a = x_before.toarray() if scipy.sparse.issparse(x_before) else np.asarray(x_before)
    b = adata_out.X.toarray() if scipy.sparse.issparse(adata_out.X) else np.asarray(adata_out.X)
    assert np.allclose(a, b)


def test_celltypist_preserves_condition() -> None:
    model = _load_model_or_skip()
    adata = _clustered(_log_normalized_with_model_genes(model))
    before = list(adata.obs["condition"])
    out, _, _ = run_primitive("celltypist", adata, model=model, majority_voting=False)
    assert list(out.obs["condition"]) == before


def test_celltypist_requires_log_normalized(adata: AnnData) -> None:
    """Raw-counts state is rejected by the requires.layer_state contract.

    This needs no model download, so it runs even offline.
    """
    from lattice_primitives.registry import SchemaViolation

    adata.uns["_lattice_layer_state"] = "raw_counts"
    with pytest.raises(SchemaViolation, match="layer_state"):
        run_primitive("celltypist", adata, model=_MODEL_NAME)

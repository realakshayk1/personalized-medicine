"""End-to-end test of the full M1 clustering chain via run_primitive.

Exercises the canonical dim_reduce -> cluster -> plot sequence on log-normalized
data, asserting state flows correctly and the condition column survives intact.
"""

from __future__ import annotations

from anndata import AnnData

import lattice_primitives.all_primitives  # noqa: F401  (register all)
from lattice_primitives import run_primitive

CHAIN = [
    ("highly_variable_genes", {"n_top_genes": 200}),
    ("pca", {"n_comps": 30, "random_state": 0}),
    ("neighbors", {"n_neighbors": 15, "random_state": 0}),
    ("umap", {"random_state": 0}),
    ("leiden", {"resolution": 1.0, "random_state": 0}),
    ("plot_umap", {"color": "leiden"}),
]


def test_full_chain_runs(adata_log_normalized: AnnData) -> None:
    adata = adata_log_normalized
    condition_before = list(adata.obs["condition"])
    figures_seen = 0

    for name, params in CHAIN:
        adata, warnings, result = run_primitive(name, adata, **params)
        # No error-severity warnings should fire on clean structured data.
        assert not [w for w in warnings if w.severity == "error"], (
            f"{name} produced error warnings: {warnings}"
        )
        figures_seen += len(result.figures)

    # Final state has the full embedding + clustering present.
    assert "highly_variable" in adata.var.columns
    assert "X_pca" in adata.obsm
    assert "X_umap" in adata.obsm
    assert "leiden" in adata.obs.columns
    assert adata.obs["leiden"].nunique() >= 2
    # One figure from plot_umap.
    assert figures_seen == 1
    # case/control protection held across the whole chain.
    assert list(adata.obs["condition"]) == condition_before

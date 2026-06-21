"""Tests for the pseudobulk aggregation primitive.

Local fixtures only (per TRACK-M6 isolation rules — nothing added to conftest).
Exercises the real run_primitive path.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from anndata import AnnData

import lattice_primitives.preprocess.pseudobulk  # noqa: F401  (register)
from lattice_primitives import run_primitive

N_SAMPLES = 6
N_GROUPS = 1
CELLS_PER_SAMPLE = 40
N_GENES = 80


def _make_raw_multisample(
    seed: int = 0,
    cells_per_sample: int = CELLS_PER_SAMPLE,
) -> AnnData:
    """Raw-count AnnData: 6 samples (3 treated / 3 control), one cell-type group.

    A block of marker genes is up-regulated in treated samples so downstream DE
    has real signal. layer_state='raw_counts'. obs has 'sample', 'condition',
    'cell_type'.
    """
    rng = np.random.default_rng(seed)
    samples = [f"s{i}" for i in range(1, N_SAMPLES + 1)]
    conds = ["treated"] * (N_SAMPLES // 2) + ["control"] * (N_SAMPLES - N_SAMPLES // 2)

    blocks = []
    obs_rows = []
    for s, c in zip(samples, conds):
        base = rng.poisson(lam=5, size=(cells_per_sample, N_GENES)).astype(np.float32)
        if c == "treated":
            base[:, :10] += rng.poisson(lam=20, size=(cells_per_sample, 10)).astype(np.float32)
        blocks.append(base)
        obs_rows.extend(
            {"sample": s, "condition": c, "cell_type": "Tcell"} for _ in range(cells_per_sample)
        )

    X = scipy.sparse.csr_matrix(np.vstack(blocks))
    adata = AnnData(X=X)
    adata.obs = pd.DataFrame(obs_rows)
    adata.obs_names = [f"CELL_{i:05d}" for i in range(adata.n_obs)]
    adata.var_names = [f"GENE{i:04d}" for i in range(N_GENES)]
    adata.uns["_lattice_layer_state"] = "raw_counts"
    return adata


@pytest.fixture
def adata_raw_multisample() -> AnnData:
    return _make_raw_multisample()


def test_pseudobulk_aggregates_to_sample_by_group(adata_raw_multisample: AnnData) -> None:
    out, warnings, result = run_primitive(
        "pseudobulk",
        adata_raw_multisample,
        sample_col="sample",
        groups_col="cell_type",
        mode="sum",
        min_cells=10,
    )
    # 6 samples x 1 group, none filtered (40 cells each >= 10).
    assert out.n_obs == N_SAMPLES * N_GROUPS
    assert out.n_vars == N_GENES
    assert out.uns["_lattice_layer_state"] == "pseudobulk_counts"
    assert result.primitive == "pseudobulk"
    # condition carried through aggregation.
    assert "condition" in out.obs.columns
    assert set(out.obs["condition"].astype(str)) == {"treated", "control"}


def test_pseudobulk_counts_are_integer_sums(adata_raw_multisample: AnnData) -> None:
    out, _, _ = run_primitive(
        "pseudobulk", adata_raw_multisample, sample_col="sample", groups_col="cell_type"
    )
    X = out.X
    x_arr = np.asarray(X.todense()) if scipy.sparse.issparse(X) else np.asarray(X)
    # Summed raw counts must be integer-valued.
    assert np.allclose(x_arr, np.round(x_arr))
    assert x_arr.min() >= 0


def test_pseudobulk_min_cells_filters_small_samples() -> None:
    # One sample with only 5 cells should be dropped at min_cells=10.
    adata = _make_raw_multisample()
    # Shrink sample 's1' to 5 cells by dropping rows.
    s1_idx = np.where(adata.obs["sample"].to_numpy() == "s1")[0]
    drop = s1_idx[5:]  # keep only first 5 of s1
    keep_mask = np.ones(adata.n_obs, dtype=bool)
    keep_mask[drop] = False
    adata = adata[keep_mask].copy()
    adata.uns["_lattice_layer_state"] = "raw_counts"

    out, _, _ = run_primitive(
        "pseudobulk", adata, sample_col="sample", groups_col="cell_type", min_cells=10
    )
    # s1 (5 cells) dropped -> 5 remaining samples.
    assert out.n_obs == N_SAMPLES - 1
    assert "s1_Tcell" not in list(out.obs_names)


def test_pseudobulk_missing_sample_col_raises(adata_raw_multisample: AnnData) -> None:
    with pytest.raises(ValueError, match="sample_col 'nope'"):
        run_primitive("pseudobulk", adata_raw_multisample, sample_col="nope")


def test_pseudobulk_missing_groups_col_raises(adata_raw_multisample: AnnData) -> None:
    with pytest.raises(ValueError, match="groups_col 'nope'"):
        run_primitive("pseudobulk", adata_raw_multisample, sample_col="sample", groups_col="nope")


def test_pseudobulk_requires_raw_counts_state(adata_raw_multisample: AnnData) -> None:
    from lattice_primitives.registry import SchemaViolation

    adata_raw_multisample.uns["_lattice_layer_state"] = "log_normalized"
    with pytest.raises(SchemaViolation, match="layer_state"):
        run_primitive(
            "pseudobulk", adata_raw_multisample, sample_col="sample", groups_col="cell_type"
        )


def test_pseudobulk_sanity_warns_on_too_few_samples() -> None:
    # Only 2 samples -> "< 4 pseudobulk samples" warning.
    rng = np.random.default_rng(1)
    blocks, obs_rows = [], []
    for s, c in [("s1", "treated"), ("s2", "control")]:
        blocks.append(rng.poisson(lam=5, size=(30, N_GENES)).astype(np.float32))
        obs_rows.extend({"sample": s, "condition": c, "cell_type": "Tcell"} for _ in range(30))
    adata = AnnData(X=scipy.sparse.csr_matrix(np.vstack(blocks)))
    adata.obs = pd.DataFrame(obs_rows)
    adata.obs_names = [f"C{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"GENE{i:04d}" for i in range(N_GENES)]
    adata.uns["_lattice_layer_state"] = "raw_counts"

    out, warnings, _ = run_primitive(
        "pseudobulk", adata, sample_col="sample", groups_col="cell_type"
    )
    assert out.n_obs == 2
    assert any("pseudobulk sample" in w.message for w in warnings)

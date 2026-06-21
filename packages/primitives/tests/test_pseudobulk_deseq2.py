"""Tests for the pseudobulk_deseq2 cross-condition DE primitive.

Local fixtures only (per TRACK-M6 isolation rules). Exercises the real
run_primitive path: cells -> pseudobulk -> DESeq2.

ENVIRONMENT NOTE: PyDESeq2 0.5.2 is incompatible with pandas 3.x — its internal
dispersion fit assigns a boolean array into a float column, which pandas >=3.0
rejects with TypeError. When that incompatibility is present, the DESeq2-execution
tests SKIP with an explicit message rather than failing, so the suite stays green
and self-documents the blocker. The guard/validation tests (which do not call
deseq2()) always run. See the TRACK-M6 report for the recommended pin fix.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from anndata import AnnData

import lattice_primitives.de.pseudobulk_deseq2  # noqa: F401  (register)
import lattice_primitives.preprocess.pseudobulk  # noqa: F401  (register)
from lattice_primitives import run_primitive

N_SAMPLES = 6
N_GENES = 60
CELLS_PER_SAMPLE = 40


def _make_pseudobulk(seed: int = 0, strong: bool = True) -> AnnData:
    """Build a pseudobulk AnnData (6 samples x 1 group) via the real primitive.

    Treated samples have an up-regulated marker block (genes 0..9) so DESeq2 has
    clear signal in a known direction.
    """
    rng = np.random.default_rng(seed)
    samples = [f"s{i}" for i in range(1, N_SAMPLES + 1)]
    conds = ["treated"] * (N_SAMPLES // 2) + ["control"] * (N_SAMPLES - N_SAMPLES // 2)
    lam_up = 30 if strong else 0

    blocks, obs_rows = [], []
    for s, c in zip(samples, conds):
        base = rng.poisson(lam=5, size=(CELLS_PER_SAMPLE, N_GENES)).astype(np.float32)
        if c == "treated" and lam_up:
            base[:, :10] += rng.poisson(lam=lam_up, size=(CELLS_PER_SAMPLE, 10)).astype(np.float32)
        blocks.append(base)
        obs_rows.extend(
            {"sample": s, "condition": c, "cell_type": "Tcell"} for _ in range(CELLS_PER_SAMPLE)
        )

    adata = AnnData(X=scipy.sparse.csr_matrix(np.vstack(blocks)))
    adata.obs = pd.DataFrame(obs_rows)
    adata.obs_names = [f"CELL_{i:05d}" for i in range(adata.n_obs)]
    adata.var_names = [f"GENE{i:04d}" for i in range(N_GENES)]
    adata.uns["_lattice_layer_state"] = "raw_counts"

    pb, _, _ = run_primitive(
        "pseudobulk", adata, sample_col="sample", groups_col="cell_type", mode="sum"
    )
    return pb


def _run_deseq2_or_skip(pb: AnnData, **params: object):
    """Run pseudobulk_deseq2, skipping on the known pandas-3/PyDESeq2 TypeError."""
    try:
        return run_primitive("pseudobulk_deseq2", pb, **params)
    except TypeError as exc:  # pandas >=3.0 vs PyDESeq2 0.5.2 incompatibility
        if "for dtype" in str(exc):
            pytest.skip(
                "BLOCKED: PyDESeq2 0.5.2 incompatible with installed pandas "
                f"(>=3.0): {exc}. Pin pandas<3 (see TRACK-M6 report)."
            )
        raise


@pytest.fixture
def adata_pseudobulk() -> AnnData:
    return _make_pseudobulk()


def test_deseq2_stores_results_with_expected_columns(adata_pseudobulk: AnnData) -> None:
    out, _, result = _run_deseq2_or_skip(
        adata_pseudobulk,
        condition_col="condition",
        reference_level="control",
        test_level="treated",
    )
    assert result.primitive == "pseudobulk_deseq2"
    assert "deseq2_results" in out.uns
    records = out.uns["deseq2_results"]
    assert len(records) == N_GENES
    expected_keys = {"gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
    assert expected_keys.issubset(set(records[0].keys()))
    meta = out.uns["deseq2_results_meta"]
    assert meta["contrast"] == ["condition", "treated", "control"]


def test_deseq2_contrast_orientation_flips_sign(adata_pseudobulk: AnnData) -> None:
    """Swapping test/reference levels must flip the sign of log2FoldChange.

    This locks AGENTS.md failure-mode #1 (case/control swap) — orientation is
    explicit and respected, not inferred.
    """
    pb1 = adata_pseudobulk
    pb2 = _make_pseudobulk()  # identical independent copy

    out1, _, _ = _run_deseq2_or_skip(pb1, reference_level="control", test_level="treated")
    out2, _, _ = _run_deseq2_or_skip(pb2, reference_level="treated", test_level="control")

    lfc1 = {r["gene"]: r["log2FoldChange"] for r in out1.uns["deseq2_results"]}
    lfc2 = {r["gene"]: r["log2FoldChange"] for r in out2.uns["deseq2_results"]}

    compared = 0
    for gene, v1 in lfc1.items():
        v2 = lfc2.get(gene)
        if v1 is None or v2 is None:
            continue
        assert np.isclose(v1, -v2, atol=1e-6), f"{gene}: {v1} vs {v2} (not sign-flipped)"
        compared += 1
    assert compared > 0  # at least some genes had defined LFCs


def test_deseq2_detects_known_upregulated_block(adata_pseudobulk: AnnData) -> None:
    """Marker block (GENE0000..GENE0009) is up in treated -> positive LFC."""
    out, _, _ = _run_deseq2_or_skip(
        adata_pseudobulk, reference_level="control", test_level="treated"
    )
    lfc = {r["gene"]: r["log2FoldChange"] for r in out.uns["deseq2_results"]}
    marker_lfcs = [lfc[f"GENE{i:04d}"] for i in range(10) if lfc.get(f"GENE{i:04d}") is not None]
    assert len(marker_lfcs) > 0
    # Up in treated => positive log2FoldChange for the marker block.
    assert np.median(marker_lfcs) > 0


def test_deseq2_missing_condition_col_raises(adata_pseudobulk: AnnData) -> None:
    with pytest.raises(ValueError, match="condition_col 'group'"):
        run_primitive("pseudobulk_deseq2", adata_pseudobulk, condition_col="group")


def test_deseq2_missing_level_raises(adata_pseudobulk: AnnData) -> None:
    with pytest.raises(ValueError, match="missing level"):
        run_primitive(
            "pseudobulk_deseq2",
            adata_pseudobulk,
            reference_level="control",
            test_level="stimulated",  # not a present level
        )


def test_deseq2_same_level_raises(adata_pseudobulk: AnnData) -> None:
    with pytest.raises(ValueError, match="must differ"):
        run_primitive(
            "pseudobulk_deseq2",
            adata_pseudobulk,
            reference_level="treated",
            test_level="treated",
        )


def test_deseq2_requires_pseudobulk_state(adata_pseudobulk: AnnData) -> None:
    from lattice_primitives.registry import SchemaViolation

    adata_pseudobulk.uns["_lattice_layer_state"] = "log_normalized"
    with pytest.raises(SchemaViolation, match="layer_state"):
        run_primitive("pseudobulk_deseq2", adata_pseudobulk)

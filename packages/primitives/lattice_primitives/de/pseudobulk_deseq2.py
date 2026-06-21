"""Primitive: pseudobulk_deseq2.

Wraps PyDESeq2 (the Python port of DESeq2) to run cross-condition differential
expression on PSEUDOBULK profiles. This is the cross-condition DE path:
``pseudobulk`` (aggregate cells -> sample x group) then ``pseudobulk_deseq2``
(negative-binomial GLM with DESeq2's size-factor normalization and dispersion
shrinkage). It is distinct from the single-cell one-vs-rest marker test in
``rank_genes_groups``.

Requires layer_state="pseudobulk_counts" (produced by the ``pseudobulk``
primitive). ``condition_col`` is a PARAMETER, so it is guarded in the body, not
in static ``requires.obs_columns``.

CRITICAL — case/control orientation (AGENTS.md failure mode #1):
    The DESeq2 contrast is ALWAYS specified explicitly and oriented as
        contrast = [condition_col, test_level, reference_level]
    so a positive log2FoldChange means "up in `test_level` relative to
    `reference_level`". SWAPPING test_level and reference_level FLIPS the sign of
    every log2FoldChange and inverts which genes read as up/down-regulated. The
    orientation is therefore made explicit and never inferred. A test asserts the
    sign flips when the levels are swapped, locking this contract in place.

Additions to AnnData:
    .uns["deseq2_results"] — the DESeq2 results table as a list of per-gene
        dict records (JSON-serializable): keys
        {gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj}, plus
        .uns["deseq2_results_meta"] recording the oriented contrast for provenance.

NOTE (biology-semantics, provisional — FLAGGED FOR LSM-FOUNDER REVIEW):
    - design="~condition" (single-factor model) is the simplest valid design.
      Real studies often need to control for batch/donor (e.g. "~batch + condition");
      provisional default, flagged for LSM review.
    - reference_level="control" / test_level="treated" are provisional naming
      conventions matched to the fixture/condition vocabulary; flagged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import AnnDataState, ParamSpec, Warning_
from lattice_primitives.registry import primitive, sanity_check

# ---------------------------------------------------------------------------
# Primitive registration
# ---------------------------------------------------------------------------

_REQUIRES = AnnDataState(
    layer_state="pseudobulk_counts",
    # condition_col is a parameter; guarded in the body.
)
_PRODUCES = AnnDataState(
    adds_uns=["deseq2_results"],
    immutable_obs_columns=["condition"],  # case/control protection
)
_PARAMS = {
    "condition_col": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "condition",
            "description": (
                "The .obs column holding the experimental condition. Must be an "
                "existing column with at least two levels including both "
                "reference_level and test_level."
            ),
        }
    ),
    "reference_level": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "control",
            "description": (
                "The baseline level of condition_col (the denominator of the "
                "contrast). A positive log2FoldChange is UP relative to this."
            ),
        }
    ),
    "test_level": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "treated",
            "description": (
                "The level of condition_col being tested (the numerator of the "
                "contrast). A positive log2FoldChange is UP in this level."
            ),
        }
    ),
    "design": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "~condition",
            "description": (
                "DESeq2 design formula. Defaults to a single-factor model on "
                "condition_col. More complex designs (e.g. '~batch + condition') "
                "must include condition_col as the last term."
            ),
        }
    ),
}

# Columns of the PyDESeq2 results_df we surface (in this order).
_RESULT_COLUMNS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


@primitive(
    name="pseudobulk_deseq2",
    category="de",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",  # scanpy unused here; pin tracked for consistency
    citation_key="love_deseq2_2014",
    row_modifying=False,
)
def pseudobulk_deseq2(
    adata: AnnData,
    condition_col: str = "condition",
    reference_level: str = "control",
    test_level: str = "treated",
    design: str = "~condition",
) -> AnnData:
    """Run PyDESeq2 differential expression on pseudobulk counts.

    Parameters
    ----------
    adata:
        AnnData in pseudobulk_counts layer_state (rows are sample x group
        profiles, .X is summed raw integer counts).
    condition_col:
        .obs column holding the experimental condition.
    reference_level:
        Baseline (denominator) level of the contrast.
    test_level:
        Tested (numerator) level of the contrast.
    design:
        DESeq2 design formula string.

    Returns
    -------
    AnnData with .uns['deseq2_results'] (list of per-gene dict records) and
    .uns['deseq2_results_meta'] (the oriented contrast) added in-place.
    """
    # Imported lazily so importing the primitive module is cheap and does not
    # pull PyDESeq2's heavy dependency tree unless the primitive actually runs.
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if condition_col not in adata.obs.columns:
        raise ValueError(
            f"pseudobulk_deseq2 requires condition_col '{condition_col}' in "
            f"adata.obs. Available columns: {list(adata.obs.columns)}."
        )

    levels = set(adata.obs[condition_col].astype(str).unique())
    missing = {reference_level, test_level} - levels
    if missing:
        raise ValueError(
            f"pseudobulk_deseq2 condition_col '{condition_col}' is missing "
            f"level(s) {sorted(missing)}. Present levels: {sorted(levels)}. "
            "Both reference_level and test_level must be present."
        )
    if reference_level == test_level:
        raise ValueError(
            "pseudobulk_deseq2 reference_level and test_level must differ "
            f"(both are '{reference_level}')."
        )

    # Build an integer counts frame (samples x genes) for PyDESeq2.
    X = adata.X
    if scipy.sparse.issparse(X):
        x_arr = np.asarray(X.todense())
    else:
        x_arr = np.asarray(X)
    counts = pd.DataFrame(
        np.rint(x_arr).astype(np.int64),
        index=adata.obs_names.astype(str),
        columns=adata.var_names.astype(str),
    )
    metadata = pd.DataFrame(
        {condition_col: adata.obs[condition_col].astype(str).to_numpy()},
        index=adata.obs_names.astype(str),
    )

    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        quiet=True,
    )
    dds.deseq2()

    # CRITICAL: contrast orientation is explicit. [factor, numerator, denominator]
    # => positive log2FoldChange == up in test_level vs reference_level.
    contrast = [condition_col, test_level, reference_level]
    stats = DeseqStats(dds, contrast=contrast, quiet=True)
    stats.summary()
    results_df = stats.results_df

    # Serialize results as a list of plain dict records (JSON-safe).
    records: list[dict[str, object]] = []
    for gene, row in results_df.iterrows():
        rec: dict[str, object] = {"gene": str(gene)}
        for col in _RESULT_COLUMNS:
            val = row.get(col)
            rec[col] = (
                None if val is None or (isinstance(val, float) and np.isnan(val)) else float(val)
            )
        records.append(rec)

    adata.uns["deseq2_results"] = records
    adata.uns["deseq2_results_meta"] = {
        "condition_col": condition_col,
        "reference_level": reference_level,
        "test_level": test_level,
        "design": design,
        "contrast": contrast,
        "n_genes": len(records),
    }

    n_sig = sum(
        1
        for r in records
        if r["padj"] is not None and float(r["padj"]) < 0.05  # type: ignore[arg-type]
    )
    logger.info(
        f"pseudobulk_deseq2: design='{design}', contrast={contrast} -> "
        f"{len(records)} genes tested, {n_sig} significant (padj<0.05)"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("pseudobulk_deseq2")
def _check_pseudobulk_deseq2(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate pseudobulk DESeq2 output.

    Checks:
    1. layer_state is 'pseudobulk_counts' (DE was run on aggregated counts).
    2. Both contrast levels have >= 2 replicates (DESeq2 dispersion needs them).
    3. Results were produced and aren't entirely NaN-padj (degenerate fit).
    """
    warnings: list[Warning_] = []
    condition_col = params.get("condition_col", "condition")
    reference_level = params.get("reference_level", "control")
    test_level = params.get("test_level", "treated")

    # 1: ran on pseudobulk
    layer_state = adata_after.uns.get("_lattice_layer_state")
    if layer_state != "pseudobulk_counts":
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"layer_state is '{layer_state}', expected 'pseudobulk_counts'. "
                    "Pseudobulk DESeq2 must run on aggregated sample-level counts; "
                    "running it on per-cell data treats cells as replicates and "
                    "drastically inflates false positives."
                ),
                primitive="pseudobulk_deseq2",
            )
        )

    # 2: replication on both sides of the contrast
    if condition_col in adata_after.obs.columns:
        counts = adata_after.obs[condition_col].astype(str).value_counts()
        for level in (reference_level, test_level):
            n = int(counts.get(level, 0))
            if n < 2:
                warnings.append(
                    Warning_(
                        severity="warn",
                        message=(
                            f"Condition level '{level}' has only {n} pseudobulk "
                            "replicate(s). DESeq2 dispersion estimation is "
                            "unreliable with fewer than 2 replicates per group."
                        ),
                        primitive="pseudobulk_deseq2",
                    )
                )

    # 3: results present and not entirely NaN
    results = adata_after.uns.get("deseq2_results")
    if not results:
        warnings.append(
            Warning_(
                severity="error",
                message=(
                    "DESeq2 produced no results table. The fit may have failed; "
                    "check that counts are integer and the design is valid."
                ),
                primitive="pseudobulk_deseq2",
            )
        )
    else:
        padj_vals = [r.get("padj") for r in results]
        if all(p is None for p in padj_vals):
            warnings.append(
                Warning_(
                    severity="error",
                    message=(
                        "All adjusted p-values are NaN. DESeq2 produced no usable "
                        "results — check for zero-count genes or insufficient "
                        "replication."
                    ),
                    primitive="pseudobulk_deseq2",
                )
            )

    return warnings

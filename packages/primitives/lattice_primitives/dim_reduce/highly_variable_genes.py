"""Primitive: highly_variable_genes.

Wraps sc.pp.highly_variable_genes to flag the most variable genes for downstream
dimensionality reduction. Feature selection focuses PCA on the genes that carry
biological signal and suppresses technical noise from uniformly-expressed genes.

Default flavor is "seurat", which operates on the LOG-NORMALIZED matrix (the state
produced by normalize_total_log1p) — so this primitive runs after normalization.

NOTE (biology-semantics, provisional — see docs/BIOLOGY_DEFAULTS.md):
    flavor="seurat" + n_top_genes=2000 are provisional defaults pending LSM-founder
    review. The main alternative, flavor="seurat_v3", requires RAW counts (and would
    have to run before normalization) plus scikit-misc.

Columns added to .var (seurat flavor):
    highly_variable      — bool, True for selected genes
    means                — mean expression per gene
    dispersions          — dispersion per gene
    dispersions_norm     — normalized dispersion (selection criterion)

Layer state is passthrough (log_normalized in, log_normalized out).
"""

from __future__ import annotations

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
)
_PRODUCES = AnnDataState(
    var_columns=["highly_variable", "means", "dispersions", "dispersions_norm"],
    immutable_obs_columns=["condition"],
)
_PARAMS = {
    "n_top_genes": ParamSpec.model_validate(
        {
            "type": "int",
            "default": 2000,
            "description": (
                "Number of highly variable genes to select. 2000 is the field "
                "standard. Selection is by normalized dispersion (seurat flavor)."
            ),
            "min_value": 1,
        }
    ),
    "flavor": ParamSpec.model_validate(
        {
            "type": "str",
            "default": "seurat",
            "description": (
                "HVG method. 'seurat' (default) operates on log-normalized data. "
                "'cell_ranger' is a binned-dispersion alternative. 'seurat_v3' "
                "requires raw counts + scikit-misc and is not wired in v0.1."
            ),
            "choices": ["seurat", "cell_ranger"],
        }
    ),
    "batch_key": ParamSpec.model_validate(
        {
            "type": "str",
            "default": None,
            "description": (
                "Optional .obs column. If set, HVGs are computed per batch then "
                "combined, reducing batch-driven feature selection."
            ),
        }
    ),
}


@primitive(
    name="highly_variable_genes",
    category="dim_reduce",
    requires=_REQUIRES,
    produces=_PRODUCES,
    params=_PARAMS,
    scanpy_version=">=1.11.0,<1.12",
    citation_key="satija_hvg_2015",
    row_modifying=False,
)
def highly_variable_genes(
    adata: AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
    batch_key: str | None = None,
) -> AnnData:
    """Flag highly variable genes via sc.pp.highly_variable_genes.

    Operates on the log-normalized matrix (seurat / cell_ranger flavors). Adds
    boolean .var['highly_variable'] plus the dispersion statistics used for
    selection. Deterministic — no random_state.

    Parameters
    ----------
    adata:
        AnnData in log_normalized layer_state.
    n_top_genes:
        Number of genes to flag as highly variable.
    flavor:
        'seurat' or 'cell_ranger' (both operate on log data).
    batch_key:
        Optional .obs column for per-batch HVG selection.

    Returns
    -------
    AnnData with .var HVG columns added in-place.
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        batch_key=batch_key,
        subset=False,  # never drop genes here; PCA masks via highly_variable
        inplace=True,
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    logger.info(
        f"highly_variable_genes: flavor={flavor}, n_top_genes={n_top_genes} -> "
        f"{n_hvg} genes flagged highly_variable"
    )
    return adata


# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------


@sanity_check("highly_variable_genes")
def _check_hvg(
    adata_before: AnnData,
    adata_after: AnnData,
    params: dict,  # type: ignore[type-arg]
) -> list[Warning_]:
    """Validate HVG selection.

    Checks:
    1. Requested n_top_genes exceeds the number of genes available.
    2. Far fewer genes flagged than requested (degenerate / low-complexity input).
    """
    warnings: list[Warning_] = []

    n_requested = int(params.get("n_top_genes", 2000))
    n_genes = adata_after.n_vars
    n_flagged = int(adata_after.var["highly_variable"].sum())

    if n_requested > n_genes:
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Requested n_top_genes={n_requested} but the dataset has only "
                    f"{n_genes} genes. All genes were flagged highly_variable, so "
                    "feature selection had no effect."
                ),
                primitive="highly_variable_genes",
            )
        )
    elif n_flagged < 0.5 * min(n_requested, n_genes):
        warnings.append(
            Warning_(
                severity="warn",
                message=(
                    f"Only {n_flagged} genes flagged highly_variable (requested "
                    f"{n_requested}). This can indicate low-complexity data or that "
                    "the input was not log-normalized as the 'seurat' flavor expects."
                ),
                primitive="highly_variable_genes",
            )
        )

    return warnings

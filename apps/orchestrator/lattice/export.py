"""Standalone notebook / report export (TRACK-M4).

Given a session's provenance log, build a *runnable* Jupyter notebook that
reproduces the analysis end-to-end. The notebook is strictly primitive-based:
every analysis cell calls ``run_primitive(...)`` from the version-pinned
``lattice_primitives`` library. We never emit free-form Scanpy/Seurat code
(ADR-001 — the LLM never writes analysis code, and neither does this exporter).

Public API
----------
build_notebook(provenance_log, citations, primitive_citations) -> NotebookNode
    Pure function: turns a provenance log dict + citation metadata into an
    nbformat NotebookNode. No live server / DB required — unit-testable.

notebook_to_py(nb) -> str
    Convert a NotebookNode to a ``py:percent`` script string via jupytext.

router : APIRouter
    FastAPI router exposing ``GET /export/{session_id}?format=ipynb|py``.

The route reads the live ProvenanceStore + the primitive registry + citations
to assemble the inputs, then calls the two pure functions above.
"""

from __future__ import annotations

from typing import Any

import nbformat
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import Response
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

# ---------------------------------------------------------------------------
# Notebook builder (pure)
# ---------------------------------------------------------------------------

_LOADER_TEMPLATE = """\
import hashlib

import anndata as ad
from lattice_primitives.registry import run_primitive

# Path to the input AnnData. Place the original .h5ad next to this notebook.
INPUT_PATH = {input_path!r}

# SHA256 of the exact input used in the original Lattice session. The export is
# reproducible only against the identical input file (AGENTS.md: reproducibility
# is non-negotiable). This guard fails loudly if the file differs.
EXPECTED_SHA256 = {expected_sha256!r}

with open(INPUT_PATH, "rb") as _fh:
    _digest = hashlib.sha256(_fh.read()).hexdigest()
if EXPECTED_SHA256 and _digest != EXPECTED_SHA256:
    raise ValueError(
        f"Input hash mismatch: expected {{EXPECTED_SHA256}}, got {{_digest}}. "
        "This notebook was generated for a different input file."
    )

adata = ad.read_h5ad(INPUT_PATH)
adata.var_names_make_unique()
"""


def _format_params(params: dict[str, Any]) -> str:
    """Render primitive params as Python keyword arguments.

    Uses ``repr`` so ``None`` -> ``None``, strings get quoted, etc. — producing
    valid, copy-pasteable Python.
    """
    return ", ".join(f"{k}={v!r}" for k, v in params.items())


def _primitive_call(primitive: str, params: dict[str, Any]) -> str:
    """Reconstruct the exact primitive invocation for a code cell.

    ``run_primitive`` returns ``(adata, warnings, result)``; we rebind ``adata``
    so the pipeline threads through cell to cell.
    """
    param_str = _format_params(params)
    args = f'"{primitive}", adata'
    if param_str:
        args = f"{args}, {param_str}"
    return f"adata, _warnings, _result = run_primitive({args})"


def build_notebook(
    provenance_log: dict[str, Any],
    citations: dict[str, dict[str, Any]],
    primitive_citations: dict[str, str] | None = None,
) -> nbformat.NotebookNode:
    """Build a runnable notebook from a provenance log.

    Parameters
    ----------
    provenance_log:
        A dict mirroring ``ProvenanceLog`` (``session_id``, ``started_at``,
        ``input_data`` with ``filename``/``sha256``/``n_obs``/``n_vars``, and a
        list of ``steps`` each with ``primitive``, ``params``,
        ``user_explanation``, ``methods_paragraph``, ``input_hash``, ...).
    citations:
        Mapping of citation_key -> citation metadata (mirrors citations.yaml).
    primitive_citations:
        Mapping of primitive name -> citation_key. Used to resolve which
        citations to include in the bibliography. If omitted, the bibliography
        lists every citation key present in ``citations`` (still deterministic
        when paired with the registry-derived map).

    Cell layout
    -----------
    1 header markdown cell
    + 1 loader code cell
    + for each step: 1 markdown cell + 1 code cell  (2 * n_steps)
    + 1 aggregated "Methods" markdown cell
    + 1 bibliography markdown cell
    """
    primitive_citations = primitive_citations or {}

    input_data = provenance_log.get("input_data", {}) or {}
    steps = provenance_log.get("steps", []) or []
    session_id = provenance_log.get("session_id", "unknown")
    started_at = provenance_log.get("started_at", "")
    filename = input_data.get("filename") or "data.h5ad"
    sha256 = input_data.get("sha256") or ""
    n_obs = input_data.get("n_obs", 0)
    n_vars = input_data.get("n_vars", 0)

    cells: list[Any] = []

    # --- Header markdown -------------------------------------------------
    header_lines = [
        "# Lattice — Reproducible Analysis Notebook",
        "",
        f"- **Session:** `{session_id}`",
        f"- **Generated from session started:** {started_at}",
        f"- **Input file:** `{filename}`",
        f"- **Input dimensions:** {n_obs} cells × {n_vars} genes",
        f"- **Input SHA256:** `{sha256}`",
        "",
        "This notebook was exported by Lattice. Every analysis step calls a "
        "version-pinned primitive from `lattice_primitives` via "
        "`run_primitive(...)` — no analysis code was generated freehand "
        "(ADR-001). Re-running this notebook against the identical input file "
        "reproduces the original results bit-for-bit.",
        "",
        "## Setup",
    ]
    cells.append(new_markdown_cell("\n".join(header_lines)))

    # --- Loader code cell ------------------------------------------------
    # Prefer the step-0 recorded input_hash; fall back to the session sha256.
    input_hash = sha256
    if steps:
        step0_hash = steps[0].get("input_hash") or ""
        if step0_hash:
            input_hash = step0_hash
    loader_code = _LOADER_TEMPLATE.format(
        input_path=filename,
        expected_sha256=input_hash,
    )
    cells.append(new_code_cell(loader_code))

    # --- Per-step markdown + code ----------------------------------------
    for i, step in enumerate(steps):
        primitive = step.get("primitive", "")
        params = step.get("params", {}) or {}
        user_explanation = (step.get("user_explanation") or "").strip()
        methods_paragraph = (step.get("methods_paragraph") or "").strip()

        md_lines = [f"### Step {i + 1}: `{primitive}`", ""]
        if user_explanation:
            md_lines.append(user_explanation)
            md_lines.append("")
        if methods_paragraph:
            md_lines.append(f"_Methods:_ {methods_paragraph}")
        cells.append(new_markdown_cell("\n".join(md_lines).rstrip()))

        cells.append(new_code_cell(_primitive_call(primitive, params)))

    # --- Aggregated Methods markdown -------------------------------------
    methods_lines = ["## Methods", ""]
    any_methods = False
    for i, step in enumerate(steps):
        para = (step.get("methods_paragraph") or "").strip()
        if para:
            any_methods = True
            methods_lines.append(para)
            methods_lines.append("")
    if not any_methods:
        methods_lines.append("_No methods paragraphs were recorded for this session._")
    cells.append(new_markdown_cell("\n".join(methods_lines).rstrip()))

    # --- Bibliography markdown -------------------------------------------
    cells.append(
        new_markdown_cell(
            _build_bibliography_markdown(steps, citations, primitive_citations)
        )
    )

    nb = new_notebook(cells=cells)
    nb.metadata["lattice"] = {
        "session_id": session_id,
        "exporter": "lattice.export",
        "n_steps": len(steps),
    }
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    nb.metadata["language_info"] = {"name": "python"}
    return nb


def _cited_keys_in_order(
    steps: list[dict[str, Any]],
    primitive_citations: dict[str, str],
    citations: dict[str, dict[str, Any]],
) -> list[str]:
    """Resolve the ordered, de-duplicated list of citation keys actually used.

    Order follows first-appearance in the pipeline. If no primitive->citation
    mapping is available, fall back to every key in ``citations``.
    """
    keys: list[str] = []
    if primitive_citations:
        for step in steps:
            primitive = step.get("primitive", "")
            key = primitive_citations.get(primitive)
            if key and key in citations and key not in keys:
                keys.append(key)
    if not keys:
        keys = list(citations.keys())
    return keys


def _build_bibliography_markdown(
    steps: list[dict[str, Any]],
    citations: dict[str, dict[str, Any]],
    primitive_citations: dict[str, str],
) -> str:
    keys = _cited_keys_in_order(steps, primitive_citations, citations)
    lines = ["## References", ""]
    if not keys:
        lines.append("_No citations recorded._")
        return "\n".join(lines)
    for n, key in enumerate(keys, start=1):
        c = citations.get(key, {})
        authors = c.get("authors", "")
        year = c.get("year", "")
        title = c.get("title", key)
        journal = c.get("journal", "")
        url = c.get("url", "")
        parts = []
        if authors:
            parts.append(f"{authors}")
        if year:
            parts.append(f"({year})")
        entry = " ".join(parts)
        entry = f"{entry}. {title}." if entry else f"{title}."
        if journal:
            entry += f" *{journal}*."
        if url:
            entry += f" {url}"
        lines.append(f"{n}. {entry}")
    return "\n".join(lines)


def notebook_to_py(nb: nbformat.NotebookNode) -> str:
    """Convert a notebook to a ``py:percent`` script string via jupytext."""
    import jupytext

    return jupytext.writes(nb, fmt="py:percent")


# ---------------------------------------------------------------------------
# Registry / citations helpers (used by the route; kept out of pure builder)
# ---------------------------------------------------------------------------


def load_citations() -> dict[str, dict[str, Any]]:
    """Load citations.yaml from the primitives package."""
    import pathlib

    import yaml
    import lattice_primitives

    pkg_dir = pathlib.Path(lattice_primitives.__file__).parent
    cit_path = pkg_dir / "citations.yaml"
    with open(cit_path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return data


def build_primitive_citation_map() -> dict[str, str]:
    """Map primitive name -> citation_key from the live primitive registry."""
    from lattice_primitives._interface import PRIMITIVE_REGISTRY

    return {name: entry.citation_key for name, entry in PRIMITIVE_REGISTRY.items()}


# ---------------------------------------------------------------------------
# FastAPI route
# ---------------------------------------------------------------------------

router = APIRouter()


@router.get("/export/{session_id}")
async def export_session(
    session_id: str,
    format: str = Query("ipynb", pattern="^(ipynb|py)$"),
) -> Response:
    """Export a session's provenance as a runnable notebook or py:percent script.

    Imported lazily inside the function to avoid a hard dependency on app
    internals at module import time (keeps this module unit-testable in
    isolation).
    """
    from lattice.app import get_store

    store = get_store()
    prov = await store.get_provenance(session_id)
    if prov is None:
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")

    provenance_log = prov.model_dump()
    citations = load_citations()
    primitive_citations = build_primitive_citation_map()

    nb = build_notebook(provenance_log, citations, primitive_citations)

    if format == "py":
        body = notebook_to_py(nb)
        media_type = "text/x-python"
        ext = "py"
    else:
        body = nbformat.writes(nb)
        media_type = "application/x-ipynb+json"
        ext = "ipynb"

    filename = f"lattice_{session_id}.{ext}"
    return Response(
        content=body,
        media_type=media_type,
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"',
        },
    )

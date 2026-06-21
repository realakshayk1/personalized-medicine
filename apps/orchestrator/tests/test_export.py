"""Tests for the standalone notebook/report exporter (TRACK-M4)."""

from __future__ import annotations

from typing import Any

import nbformat
import pytest
import pytest_asyncio

# Ensure primitives are registered (needed for the route's registry map).
import lattice_primitives.all_primitives  # noqa: F401
from httpx import AsyncClient
from httpx._transports.asgi import ASGITransport

from lattice.export import (
    build_notebook,
    build_primitive_citation_map,
    load_citations,
    notebook_to_py,
)

# ---------------------------------------------------------------------------
# Fixtures: a fake provenance log mirroring provenance.py's schema
# ---------------------------------------------------------------------------

FAKE_CITATIONS: dict[str, dict[str, Any]] = {
    "theis_normalization_2019": {
        "title": "Current best practices in single-cell RNA-seq analysis: a tutorial",
        "authors": "Luecken, M.D. & Theis, F.J.",
        "year": 2019,
        "journal": "Molecular Systems Biology",
        "url": "https://doi.org/10.15252/msb.20188746",
    },
    "scanpy_calculate_qc_metrics": {
        "title": "SCANPY: large-scale single-cell gene expression data analysis",
        "authors": "Wolf, F.A., Angerer, P., & Theis, F.J.",
        "year": 2018,
        "journal": "Genome Biology",
        "url": "https://doi.org/10.1186/s13059-017-1382-0",
    },
    "unused_citation": {
        "title": "Should not appear in bibliography",
        "authors": "Nobody",
        "year": 2099,
    },
}

PRIMITIVE_CITATIONS = {
    "calculate_qc_metrics": "scanpy_calculate_qc_metrics",
    "normalize_total_log1p": "theis_normalization_2019",
}


def make_fake_provenance_log() -> dict[str, Any]:
    """Build a provenance-log dict matching ProvenanceLog.model_dump()."""
    return {
        "session_id": "sess-abc-123",
        "started_at": "2026-06-21T00:00:00Z",
        "input_data": {
            "filename": "pbmc3k.h5ad",
            "sha256": "a" * 64,
            "n_obs": 2700,
            "n_vars": 32738,
        },
        "steps": [
            {
                "step_id": 0,
                "started_at": "2026-06-21T00:00:01Z",
                "primitive": "calculate_qc_metrics",
                "params": {"mt_prefix": "MT-"},
                "input_hash": "a" * 64,
                "output_hash": "b" * 64,
                "duration_sec": 0.1,
                "warnings": [],
                "user_explanation": "Computed QC metrics per cell.",
                "methods_paragraph": "QC metrics were computed using Scanpy.",
                "n_obs_before": 2700,
                "n_obs_after": 2700,
                "n_vars_before": 32738,
                "n_vars_after": 32738,
                "figures": [],
            },
            {
                "step_id": 1,
                "started_at": "2026-06-21T00:00:02Z",
                "primitive": "normalize_total_log1p",
                "params": {"target_sum": None, "exclude_highly_expressed": False},
                "input_hash": "b" * 64,
                "output_hash": "c" * 64,
                "duration_sec": 0.2,
                "warnings": [],
                "user_explanation": "Normalized and log-transformed counts.",
                "methods_paragraph": (
                    "Counts were normalized per cell and log1p-transformed "
                    "(Luecken & Theis 2019)."
                ),
                "n_obs_before": 2700,
                "n_obs_after": 2700,
                "n_vars_before": 32738,
                "n_vars_after": 32738,
                "figures": [],
            },
        ],
    }


# ---------------------------------------------------------------------------
# Pure builder tests
# ---------------------------------------------------------------------------


class TestBuildNotebook:
    def test_cell_count(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        n_steps = len(log["steps"])
        # header + loader + 2*nsteps + methods + citations
        expected = 1 + 1 + 2 * n_steps + 1 + 1
        assert len(nb.cells) == expected

    def test_notebook_is_valid(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        # Raises nbformat.ValidationError if invalid.
        nbformat.validate(nb)

    def test_loader_cell_has_hash_check(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        loader = nb.cells[1]
        assert loader.cell_type == "code"
        assert "hashlib.sha256" in loader.source
        # step-0 input_hash should be embedded as the expected hash
        assert log["steps"][0]["input_hash"] in loader.source
        assert "read_h5ad" in loader.source
        assert "pbmc3k.h5ad" in loader.source

    def test_code_cells_use_run_primitive_with_params(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        code_sources = [c.source for c in nb.cells if c.cell_type == "code"]
        joined = "\n".join(code_sources)
        # Primitive-based, never free-form analysis code (ADR-001)
        assert 'run_primitive("calculate_qc_metrics", adata' in joined
        assert "mt_prefix=" in joined
        assert "'MT-'" in joined or '"MT-"' in joined
        assert 'run_primitive("normalize_total_log1p", adata' in joined
        assert "target_sum=None" in joined
        assert "exclude_highly_expressed=False" in joined

    def test_no_freeform_scanpy_calls(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        code_sources = [c.source for c in nb.cells if c.cell_type == "code"]
        # The loader uses anndata.read_h5ad; analysis cells must only call
        # run_primitive, never sc.pp / sc.tl / sc.pp.* directly.
        for src in code_sources:
            assert "sc.pp." not in src
            assert "sc.tl." not in src

    def test_markdown_cells_carry_explanations(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        md = "\n".join(c.source for c in nb.cells if c.cell_type == "markdown")
        assert "Computed QC metrics per cell." in md
        assert "Normalized and log-transformed counts." in md

    def test_aggregated_methods_section(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        md_sources = [c.source for c in nb.cells if c.cell_type == "markdown"]
        methods = [s for s in md_sources if s.startswith("## Methods")]
        assert len(methods) == 1
        assert "QC metrics were computed using Scanpy." in methods[0]

    def test_bibliography_contains_cited_titles_only(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        md_sources = [c.source for c in nb.cells if c.cell_type == "markdown"]
        bib = [s for s in md_sources if s.startswith("## References")]
        assert len(bib) == 1
        biblio = bib[0]
        assert FAKE_CITATIONS["scanpy_calculate_qc_metrics"]["title"] in biblio
        assert FAKE_CITATIONS["theis_normalization_2019"]["title"] in biblio
        # The unused citation must NOT appear
        assert FAKE_CITATIONS["unused_citation"]["title"] not in biblio

    def test_empty_provenance_log(self) -> None:
        log = {
            "session_id": "empty",
            "started_at": "2026-06-21T00:00:00Z",
            "input_data": {
                "filename": "x.h5ad",
                "sha256": "f" * 64,
                "n_obs": 0,
                "n_vars": 0,
            },
            "steps": [],
        }
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        # header + loader + methods + citations (no per-step cells)
        assert len(nb.cells) == 4
        nbformat.validate(nb)


class TestNotebookToPy:
    def test_py_export_runnable_looking(self) -> None:
        log = make_fake_provenance_log()
        nb = build_notebook(log, FAKE_CITATIONS, PRIMITIVE_CITATIONS)
        py = notebook_to_py(nb)
        assert isinstance(py, str)
        # jupytext py:percent markers
        assert "# %%" in py
        # primitive calls survive the conversion
        assert "run_primitive" in py
        assert "calculate_qc_metrics" in py
        assert "import hashlib" in py
        # py:percent should compile as valid Python source
        compile(py, "<exported>", "exec")


# ---------------------------------------------------------------------------
# Registry/citations helper tests
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_load_citations_returns_known_keys(self) -> None:
        cits = load_citations()
        assert "theis_normalization_2019" in cits
        assert "title" in cits["theis_normalization_2019"]

    def test_primitive_citation_map_resolves(self) -> None:
        pmap = build_primitive_citation_map()
        # normalize_total_log1p is the reference primitive
        assert pmap.get("normalize_total_log1p") == "theis_normalization_2019"


# ---------------------------------------------------------------------------
# Route tests against a seeded in-memory ProvenanceStore
# ---------------------------------------------------------------------------


@pytest_asyncio.fixture
async def client() -> AsyncClient:  # type: ignore[misc]
    """Async client for the app with the export router included.

    The app under test does not yet wire export.router (TRACK-M8 applies the
    two reported lines). We include it here so the route is exercisable in
    isolation without editing app.py.
    """
    from lattice.app import app, lifespan
    from lattice.export import router as export_router

    # Idempotent include — guard against double registration across tests.
    if not any(getattr(r, "path", None) == "/export/{session_id}" for r in app.routes):
        app.include_router(export_router)

    async with lifespan(app):
        async with AsyncClient(
            transport=ASGITransport(app=app), base_url="http://test"
        ) as ac:
            yield ac  # type: ignore[misc]


class TestExportRoute:
    @pytest.mark.asyncio
    async def test_export_ipynb(self, client: AsyncClient) -> None:
        # Seed a session + steps via the live store.
        from lattice.app import get_store

        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]

        store = get_store()
        await store.update_session_upload(
            session_id=session_id,
            filename="pbmc3k.h5ad",
            sha256="a" * 64,
            n_obs=2700,
            n_vars=32738,
        )
        await store.record_step(
            session_id=session_id,
            step_id=0,
            started_at="2026-06-21T00:00:01Z",
            entry={
                "primitive": "calculate_qc_metrics",
                "params": {"mt_prefix": "MT-"},
                "input_hash": "a" * 64,
                "output_hash": "b" * 64,
                "duration_sec": 0.1,
                "warnings": [],
                "user_explanation": "Computed QC metrics.",
                "methods_paragraph": "QC metrics computed with Scanpy.",
                "n_obs_before": 2700,
                "n_obs_after": 2700,
                "n_vars_before": 32738,
                "n_vars_after": 32738,
                "figures": [],
            },
        )

        resp = await client.get(f"/export/{session_id}?format=ipynb")
        assert resp.status_code == 200
        assert "ipynb" in resp.headers["content-type"]
        nb = nbformat.reads(resp.text, as_version=4)
        nbformat.validate(nb)
        joined = "\n".join(c.source for c in nb.cells)
        assert "calculate_qc_metrics" in joined

    @pytest.mark.asyncio
    async def test_export_py(self, client: AsyncClient) -> None:
        from lattice.app import get_store

        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]
        store = get_store()
        await store.update_session_upload(
            session_id=session_id,
            filename="pbmc3k.h5ad",
            sha256="a" * 64,
            n_obs=2700,
            n_vars=32738,
        )

        resp = await client.get(f"/export/{session_id}?format=py")
        assert resp.status_code == 200
        assert "python" in resp.headers["content-type"]
        assert "# %%" in resp.text
        assert "run_primitive" in resp.text or "read_h5ad" in resp.text

    @pytest.mark.asyncio
    async def test_export_missing_session(self, client: AsyncClient) -> None:
        resp = await client.get("/export/no-such-session?format=ipynb")
        assert resp.status_code == 404

    @pytest.mark.asyncio
    async def test_export_bad_format(self, client: AsyncClient) -> None:
        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]
        resp = await client.get(f"/export/{session_id}?format=docx")
        # FastAPI query pattern validation -> 422
        assert resp.status_code == 422

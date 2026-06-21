"""Lattice orchestrator — FastAPI application.

Routes
------
GET  /health
POST /sessions
POST /sessions/{session_id}/upload
POST /sessions/{session_id}/plan
POST /sessions/{session_id}/execute  (SSE stream — ADR-002)
GET  /sessions/{session_id}/provenance
"""

from __future__ import annotations

import datetime
import hashlib
import io  # noqa: F401
import os
import tempfile
import uuid
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from typing import Any

import anndata

# ---------------------------------------------------------------------------
# Ensure primitives are registered at startup
# ---------------------------------------------------------------------------
import lattice_primitives.all_primitives  # noqa: F401, E402
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from loguru import logger
from sse_starlette.sse import EventSourceResponse

from lattice.executor import run_plan
from lattice.models import (
    ExecuteRequest,
    Plan,
    PlanRequest,
    ProvenanceLog,
    SessionCreateResponse,
    UploadResponse,
)
from lattice import export as export_module
from lattice.planner import Planner
from lattice.provenance import ProvenanceStore
from lattice.sandbox import SandboxProtocol, make_sandbox

# ---------------------------------------------------------------------------
# App state
# ---------------------------------------------------------------------------

_DB_PATH = os.environ.get("LATTICE_DB_PATH", ":memory:")

_store: ProvenanceStore | None = None
# Per-session AnnData summaries stored in memory (production would use SQLite)
_session_summaries: dict[str, dict[str, Any]] = {}
# Per-session sandbox instances (FakeSandbox by default; see make_sandbox)
_session_sandboxes: dict[str, SandboxProtocol] = {}


def get_store() -> ProvenanceStore:
    if _store is None:
        raise RuntimeError("ProvenanceStore not initialized")
    return _store


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    global _store
    _store = ProvenanceStore(db_path=_DB_PATH)
    await _store.open()
    logger.info(f"ProvenanceStore opened at {_DB_PATH}")
    yield
    if _store is not None:
        await _store.close()
        logger.info("ProvenanceStore closed")


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------

app = FastAPI(
    title="Lattice Orchestrator",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # open for localhost dev; restrict in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Notebook/report export routes (/export/{session_id})
app.include_router(export_module.router)


# ---------------------------------------------------------------------------
# Health
# ---------------------------------------------------------------------------


@app.get("/health")
async def health() -> dict[str, bool]:
    return {"ok": True}


# ---------------------------------------------------------------------------
# Sessions
# ---------------------------------------------------------------------------


@app.post("/sessions", response_model=SessionCreateResponse, status_code=201)
async def create_session() -> SessionCreateResponse:
    session_id = str(uuid.uuid4())
    started_at = datetime.datetime.utcnow().isoformat() + "Z"
    store = get_store()
    await store.create_session(session_id=session_id, started_at=started_at)
    _session_summaries[session_id] = {}
    logger.info(f"Created session {session_id}")
    return SessionCreateResponse(session_id=session_id)


# ---------------------------------------------------------------------------
# Upload
# ---------------------------------------------------------------------------


@app.post("/sessions/{session_id}/upload", response_model=UploadResponse)
async def upload_h5ad(
    session_id: str,
    file: UploadFile = File(...),  # noqa: B008
) -> UploadResponse:
    store = get_store()
    if not await store.session_exists(session_id):
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")

    contents = await file.read()

    # Compute sha256
    sha256 = hashlib.sha256(contents).hexdigest()

    # Load AnnData — anndata.read_h5ad requires a path on disk, not BytesIO.
    try:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp.write(contents)
            tmp_path = tmp.name
        try:
            adata = anndata.read_h5ad(tmp_path)
        finally:
            os.unlink(tmp_path)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(
            status_code=422, detail=f"Failed to read h5ad file: {exc}"
        ) from exc

    # Ensure unique var_names (per AGENTS.md pitfall)
    adata.var_names_make_unique()

    filename = file.filename or "data.h5ad"
    n_obs = adata.n_obs
    n_vars = adata.n_vars
    obs_columns = list(adata.obs.columns)

    # Store in provenance DB
    await store.update_session_upload(
        session_id=session_id,
        filename=filename,
        sha256=sha256,
        n_obs=n_obs,
        n_vars=n_vars,
    )

    # Store summary and sandbox in memory
    _session_summaries[session_id] = {
        "n_obs": n_obs,
        "n_vars": n_vars,
        "obs_columns": obs_columns,
        "filename": filename,
    }
    sandbox = make_sandbox()
    await sandbox.upload_anndata(adata)
    _session_sandboxes[session_id] = sandbox

    logger.info(f"Session {session_id}: uploaded {filename} ({n_obs}x{n_vars})")
    return UploadResponse(
        session_id=session_id,
        filename=filename,
        sha256=sha256,
        n_obs=n_obs,
        n_vars=n_vars,
        obs_columns=obs_columns,
    )


# ---------------------------------------------------------------------------
# Plan
# ---------------------------------------------------------------------------


@app.post("/sessions/{session_id}/plan", response_model=Plan)
async def create_plan(
    session_id: str,
    body: PlanRequest,
) -> Plan:
    store = get_store()
    if not await store.session_exists(session_id):
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")

    summary = _session_summaries.get(session_id, {})
    planner = Planner()

    try:
        plan = await planner.plan(
            session_id=session_id,
            user_message=body.user_message,
            anndata_summary=summary,
        )
    except Exception as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    return plan


# ---------------------------------------------------------------------------
# Execute (SSE stream) — ADR-002: SSE chosen
# ---------------------------------------------------------------------------


@app.post("/sessions/{session_id}/execute")
async def execute_plan_endpoint(
    session_id: str,
    body: ExecuteRequest,
) -> EventSourceResponse:
    store = get_store()
    if not await store.session_exists(session_id):
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")

    sandbox = _session_sandboxes.get(session_id)
    if sandbox is None:
        # Create a fresh FakeSandbox for sessions without an upload (test convenience)
        sandbox = make_sandbox()
        _session_sandboxes[session_id] = sandbox

    plan = body.plan

    async def event_generator() -> AsyncIterator[dict[str, Any]]:
        async for event in run_plan(plan=plan, sandbox=sandbox, store=store):
            yield {"data": event.model_dump_json()}

    return EventSourceResponse(event_generator())


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


@app.get("/sessions/{session_id}/provenance", response_model=ProvenanceLog)
async def get_provenance(session_id: str) -> ProvenanceLog:
    store = get_store()
    prov = await store.get_provenance(session_id)
    if prov is None:
        raise HTTPException(status_code=404, detail=f"Session '{session_id}' not found")
    return prov

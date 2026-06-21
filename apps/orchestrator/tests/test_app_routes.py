"""Tests for FastAPI routes using TestClient."""

from __future__ import annotations

from unittest.mock import AsyncMock, patch

import anndata
import lattice_primitives.preprocess.filter_cells_basic  # noqa: F401
import lattice_primitives.preprocess.normalize_total_log1p  # noqa: F401

# Ensure primitives are imported
import lattice_primitives.qc.calculate_qc_metrics  # noqa: F401
import numpy as np
import pytest
import pytest_asyncio
import scipy.sparse
from httpx import AsyncClient
from httpx._transports.asgi import ASGITransport


def make_tiny_h5ad(n_obs: int = 10, n_vars: int = 20) -> bytes:
    import tempfile
    import pandas as pd
    anndata.settings.allow_write_nullable_strings = True
    rng = np.random.default_rng(0)
    X = scipy.sparse.csr_matrix(
        rng.negative_binomial(5, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    )
    adata = anndata.AnnData(X=X)
    adata.obs_names = pd.Index([f"cell_{i}" for i in range(n_obs)])
    adata.var_names = pd.Index([f"gene_{i}" for i in range(n_vars)])
    adata.uns["_lattice_layer_state"] = "raw_counts"
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
        tmp_path = f.name
    adata.write_h5ad(tmp_path)
    with open(tmp_path, "rb") as f:
        return f.read()


@pytest_asyncio.fixture
async def client() -> AsyncClient:  # type: ignore[misc]
    """Async TestClient for the FastAPI app (with lifespan)."""
    from lattice.app import app, lifespan

    async with lifespan(app):
        async with AsyncClient(
            transport=ASGITransport(app=app), base_url="http://test"
        ) as ac:
            yield ac  # type: ignore[misc]


class TestHealth:
    @pytest.mark.asyncio
    async def test_health_ok(self, client: AsyncClient) -> None:
        resp = await client.get("/health")
        assert resp.status_code == 200
        assert resp.json() == {"ok": True}


class TestSessions:
    @pytest.mark.asyncio
    async def test_create_session(self, client: AsyncClient) -> None:
        resp = await client.post("/sessions")
        assert resp.status_code == 201
        data = resp.json()
        assert "session_id" in data
        assert len(data["session_id"]) > 0

    @pytest.mark.asyncio
    async def test_create_multiple_sessions_unique_ids(self, client: AsyncClient) -> None:
        resp1 = await client.post("/sessions")
        resp2 = await client.post("/sessions")
        assert resp1.json()["session_id"] != resp2.json()["session_id"]


class TestUpload:
    @pytest.mark.asyncio
    async def test_upload_h5ad(self, client: AsyncClient) -> None:
        # Create session
        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]

        # Upload
        h5ad_bytes = make_tiny_h5ad()
        resp = await client.post(
            f"/sessions/{session_id}/upload",
            files={"file": ("test.h5ad", h5ad_bytes, "application/octet-stream")},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["session_id"] == session_id
        assert data["filename"] == "test.h5ad"
        assert "sha256" in data
        assert data["n_obs"] == 10
        assert data["n_vars"] == 20
        assert isinstance(data["obs_columns"], list)

    @pytest.mark.asyncio
    async def test_upload_missing_session(self, client: AsyncClient) -> None:
        h5ad_bytes = make_tiny_h5ad()
        resp = await client.post(
            "/sessions/nonexistent-session/upload",
            files={"file": ("test.h5ad", h5ad_bytes, "application/octet-stream")},
        )
        assert resp.status_code == 404

    @pytest.mark.asyncio
    async def test_upload_invalid_file(self, client: AsyncClient) -> None:
        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]

        resp = await client.post(
            f"/sessions/{session_id}/upload",
            files={"file": ("bad.h5ad", b"not valid h5ad content", "application/octet-stream")},
        )
        assert resp.status_code == 422


class TestPlan:
    @pytest.mark.asyncio
    async def test_plan_with_stub_llm(self, client: AsyncClient) -> None:
        """Test /plan endpoint using a stub that returns known JSON."""
        import json as _json

        stub_plan = {
            "session_id": "will-be-replaced",
            "steps": [
                {
                    "primitive": "calculate_qc_metrics",
                    "params": {"mt_prefix": "MT-"},
                    "rationale": "QC",
                },
            ],
            "rationale": "Test",
            "workflow": "standard_scrnaseq",
        }

        async def fake_chat(system: str, user: str) -> str:
            return _json.dumps(stub_plan)

        # Create session
        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]

        # Patch the Planner's default chat fn
        with patch("lattice.app.Planner") as MockPlanner:
            instance = MockPlanner.return_value
            instance.plan = AsyncMock(
                return_value=__import__(
                    "lattice.models", fromlist=["Plan"]
                ).Plan(
                    session_id=session_id,
                    steps=[
                        __import__(
                            "lattice.models", fromlist=["PrimitiveInvocation"]
                        ).PrimitiveInvocation(
                            primitive="calculate_qc_metrics",
                            params={"mt_prefix": "MT-"},
                        )
                    ],
                    rationale="Test",
                    workflow="standard_scrnaseq",
                )
            )

            resp = await client.post(
                f"/sessions/{session_id}/plan",
                json={"user_message": "Run QC"},
            )

        assert resp.status_code == 200
        data = resp.json()
        assert data["session_id"] == session_id
        assert len(data["steps"]) >= 1

    @pytest.mark.asyncio
    async def test_plan_missing_session(self, client: AsyncClient) -> None:
        resp = await client.post(
            "/sessions/no-such-session/plan",
            json={"user_message": "Do something"},
        )
        assert resp.status_code == 404


class TestProvenance:
    @pytest.mark.asyncio
    async def test_provenance_empty_session(self, client: AsyncClient) -> None:
        sess_resp = await client.post("/sessions")
        session_id = sess_resp.json()["session_id"]

        resp = await client.get(f"/sessions/{session_id}/provenance")
        assert resp.status_code == 200
        data = resp.json()
        assert data["session_id"] == session_id
        assert data["steps"] == []

    @pytest.mark.asyncio
    async def test_provenance_missing_session(self, client: AsyncClient) -> None:
        resp = await client.get("/sessions/nonexistent/provenance")
        assert resp.status_code == 404

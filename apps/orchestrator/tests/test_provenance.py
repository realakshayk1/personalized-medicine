"""Tests for ProvenanceStore and the schema-frozen hash."""

from __future__ import annotations

import pytest
from lattice.provenance import SCHEMA_HASH, ProvenanceStore

# ---------------------------------------------------------------------------
# Schema-frozen test — update this hash ONLY when you write a new ADR
# ---------------------------------------------------------------------------

EXPECTED_SCHEMA_HASH = "0c42781f382f36de334bfddb6772e759ba6bfe4a6d006bed2b7186952e0e28b0"


def test_provenance_schema_frozen() -> None:
    """The DDL hash must match the known value. Change requires an ADR."""
    assert SCHEMA_HASH == EXPECTED_SCHEMA_HASH, (
        f"ProvenanceStore DDL hash changed: {SCHEMA_HASH!r} != {EXPECTED_SCHEMA_HASH!r}. "
        "If you intentionally changed the schema, write a new ADR and update "
        "EXPECTED_SCHEMA_HASH in test_provenance.py."
    )


# ---------------------------------------------------------------------------
# Basic store tests
# ---------------------------------------------------------------------------


class TestProvenanceStore:
    @pytest.mark.asyncio
    async def test_create_and_exists(self, prov_store: ProvenanceStore) -> None:
        await prov_store.create_session("sess-1", "2026-05-22T00:00:00Z")
        assert await prov_store.session_exists("sess-1")

    @pytest.mark.asyncio
    async def test_not_exists(self, prov_store: ProvenanceStore) -> None:
        assert not await prov_store.session_exists("nonexistent")

    @pytest.mark.asyncio
    async def test_update_upload(self, prov_store: ProvenanceStore) -> None:
        await prov_store.create_session("sess-2", "2026-05-22T00:00:00Z")
        await prov_store.update_session_upload(
            session_id="sess-2",
            filename="data.h5ad",
            sha256="abc123",
            n_obs=100,
            n_vars=200,
        )
        prov = await prov_store.get_provenance("sess-2")
        assert prov is not None
        assert prov.input_data.filename == "data.h5ad"
        assert prov.input_data.sha256 == "abc123"
        assert prov.input_data.n_obs == 100
        assert prov.input_data.n_vars == 200

    @pytest.mark.asyncio
    async def test_record_and_fetch_steps(self, prov_store: ProvenanceStore) -> None:
        await prov_store.create_session("sess-3", "2026-05-22T00:00:00Z")

        step_entry = {
            "primitive": "calculate_qc_metrics",
            "params": {"mt_prefix": "MT-"},
            "input_hash": "aaa",
            "output_hash": "bbb",
            "duration_sec": 1.5,
            "warnings": [],
            "user_explanation": "Computed QC metrics.",
            "methods_paragraph": "QC metrics were computed using sc.pp.calculate_qc_metrics.",
            "n_obs_before": 100,
            "n_obs_after": 100,
            "n_vars_before": 200,
            "n_vars_after": 200,
            "figures": [],
        }

        await prov_store.record_step(
            session_id="sess-3",
            step_id=0,
            started_at="2026-05-22T00:01:00Z",
            entry=step_entry,
        )

        prov = await prov_store.get_provenance("sess-3")
        assert prov is not None
        assert len(prov.steps) == 1
        step = prov.steps[0]
        assert step.primitive == "calculate_qc_metrics"
        assert step.input_hash == "aaa"
        assert step.output_hash == "bbb"
        assert step.duration_sec == 1.5

    @pytest.mark.asyncio
    async def test_get_provenance_nonexistent(self, prov_store: ProvenanceStore) -> None:
        result = await prov_store.get_provenance("does-not-exist")
        assert result is None

    @pytest.mark.asyncio
    async def test_multiple_steps_ordered(self, prov_store: ProvenanceStore) -> None:
        await prov_store.create_session("sess-4", "2026-05-22T00:00:00Z")

        for i, prim in enumerate(["calculate_qc_metrics", "filter_cells_min_counts"]):
            await prov_store.record_step(
                session_id="sess-4",
                step_id=i,
                started_at=f"2026-05-22T00:0{i}:00Z",
                entry={
                    "primitive": prim,
                    "params": {},
                    "input_hash": f"in_{i}",
                    "output_hash": f"out_{i}",
                    "duration_sec": 0.5,
                    "warnings": [],
                    "user_explanation": "",
                    "methods_paragraph": "",
                    "n_obs_before": 100,
                    "n_obs_after": 90,
                    "n_vars_before": 200,
                    "n_vars_after": 200,
                    "figures": [],
                },
            )

        prov = await prov_store.get_provenance("sess-4")
        assert prov is not None
        assert len(prov.steps) == 2
        assert prov.steps[0].step_id == 0
        assert prov.steps[1].step_id == 1

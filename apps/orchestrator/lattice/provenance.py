"""Provenance store backed by SQLite via aiosqlite.

SCHEMA IS FROZEN — per AGENTS.md "never delete or rewrite the provenance log
schema." Any change requires a new ADR and an update to SCHEMA_HASH below.

The test test_provenance_schema_frozen.py asserts SCHEMA_HASH == known value.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

import aiosqlite
from loguru import logger

from lattice.models import InputData, ProvenanceEntry, ProvenanceLog

# ---------------------------------------------------------------------------
# Schema DDL — frozen; hash guards against accidental mutation
# ---------------------------------------------------------------------------

_SESSIONS_DDL = """\
CREATE TABLE IF NOT EXISTS sessions (
    session_id TEXT PRIMARY KEY,
    started_at TEXT NOT NULL,
    filename TEXT,
    sha256 TEXT,
    n_obs INTEGER,
    n_vars INTEGER
)"""

_STEPS_DDL = """\
CREATE TABLE IF NOT EXISTS steps (
    session_id TEXT NOT NULL,
    step_id INTEGER NOT NULL,
    started_at TEXT NOT NULL,
    primitive TEXT NOT NULL,
    params_json TEXT,
    input_hash TEXT,
    output_hash TEXT,
    duration_sec REAL,
    warnings_json TEXT,
    user_explanation TEXT,
    methods_paragraph TEXT,
    n_obs_before INTEGER,
    n_obs_after INTEGER,
    n_vars_before INTEGER,
    n_vars_after INTEGER,
    figures_json TEXT,
    PRIMARY KEY (session_id, step_id)
)"""

_DDL_FOR_HASH = _SESSIONS_DDL + "\n---\n" + _STEPS_DDL

SCHEMA_HASH: str = hashlib.sha256(_DDL_FOR_HASH.encode()).hexdigest()


# ---------------------------------------------------------------------------
# ProvenanceStore
# ---------------------------------------------------------------------------


class ProvenanceStore:
    """Async SQLite-backed provenance store."""

    def __init__(self, db_path: str = ":memory:") -> None:
        self._db_path = db_path
        self._db: aiosqlite.Connection | None = None

    async def open(self) -> None:
        """Open the database and create tables if needed."""
        self._db = await aiosqlite.connect(self._db_path)
        await self._db.execute("PRAGMA journal_mode=WAL")
        await self._db.execute(_SESSIONS_DDL)
        await self._db.execute(_STEPS_DDL)
        await self._db.commit()
        logger.debug(f"ProvenanceStore opened at {self._db_path}")

    async def close(self) -> None:
        if self._db is not None:
            await self._db.close()
            self._db = None

    def _conn(self) -> aiosqlite.Connection:
        if self._db is None:
            raise RuntimeError("ProvenanceStore not open — call open() first")
        return self._db

    # ------------------------------------------------------------------
    # Sessions
    # ------------------------------------------------------------------

    async def create_session(
        self,
        session_id: str,
        started_at: str,
        filename: str | None = None,
        sha256: str | None = None,
        n_obs: int | None = None,
        n_vars: int | None = None,
    ) -> None:
        await self._conn().execute(
            "INSERT OR IGNORE INTO sessions (session_id, started_at, filename, sha256, n_obs, n_vars) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (session_id, started_at, filename, sha256, n_obs, n_vars),
        )
        await self._conn().commit()

    async def update_session_upload(
        self,
        session_id: str,
        filename: str,
        sha256: str,
        n_obs: int,
        n_vars: int,
    ) -> None:
        await self._conn().execute(
            "UPDATE sessions SET filename=?, sha256=?, n_obs=?, n_vars=? WHERE session_id=?",
            (filename, sha256, n_obs, n_vars, session_id),
        )
        await self._conn().commit()

    # ------------------------------------------------------------------
    # Steps
    # ------------------------------------------------------------------

    async def record_step(
        self,
        session_id: str,
        step_id: int,
        started_at: str,
        entry: dict[str, Any],
    ) -> None:
        await self._conn().execute(
            """
            INSERT INTO steps (
                session_id, step_id, started_at, primitive, params_json,
                input_hash, output_hash, duration_sec, warnings_json,
                user_explanation, methods_paragraph,
                n_obs_before, n_obs_after, n_vars_before, n_vars_after,
                figures_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                session_id,
                step_id,
                started_at,
                entry.get("primitive", ""),
                json.dumps(entry.get("params", {})),
                entry.get("input_hash", ""),
                entry.get("output_hash", ""),
                entry.get("duration_sec", 0.0),
                json.dumps([w if isinstance(w, dict) else w for w in entry.get("warnings", [])]),
                entry.get("user_explanation", ""),
                entry.get("methods_paragraph", ""),
                entry.get("n_obs_before"),
                entry.get("n_obs_after"),
                entry.get("n_vars_before"),
                entry.get("n_vars_after"),
                json.dumps(entry.get("figures", [])),
            ),
        )
        await self._conn().commit()

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    async def get_provenance(self, session_id: str) -> ProvenanceLog | None:
        conn = self._conn()

        # Load session row
        async with conn.execute(
            "SELECT session_id, started_at, filename, sha256, n_obs, n_vars FROM sessions WHERE session_id=?",
            (session_id,),
        ) as cur:
            row = await cur.fetchone()
        if row is None:
            return None

        sess_id, started_at, filename, sha256, n_obs, n_vars = row

        # Load steps
        async with conn.execute(
            """
            SELECT step_id, started_at, primitive, params_json,
                   input_hash, output_hash, duration_sec, warnings_json,
                   user_explanation, methods_paragraph,
                   n_obs_before, n_obs_after, n_vars_before, n_vars_after,
                   figures_json
            FROM steps WHERE session_id=? ORDER BY step_id
            """,
            (session_id,),
        ) as cur:
            step_rows = await cur.fetchall()

        steps: list[ProvenanceEntry] = []
        for sr in step_rows:
            (
                step_id, s_started_at, primitive, params_json,
                input_hash, output_hash, duration_sec, warnings_json,
                user_explanation, methods_paragraph,
                n_obs_before, n_obs_after, n_vars_before, n_vars_after,
                figures_json,
            ) = sr
            warnings_raw = json.loads(warnings_json) if warnings_json else []
            steps.append(
                ProvenanceEntry(
                    step_id=step_id,
                    started_at=s_started_at,
                    primitive=primitive,
                    params=json.loads(params_json) if params_json else {},
                    input_hash=input_hash or "",
                    output_hash=output_hash or "",
                    duration_sec=duration_sec or 0.0,
                    warnings=warnings_raw,
                    user_explanation=user_explanation or "",
                    methods_paragraph=methods_paragraph or "",
                    n_obs_before=n_obs_before,
                    n_obs_after=n_obs_after,
                    n_vars_before=n_vars_before,
                    n_vars_after=n_vars_after,
                    figures=json.loads(figures_json) if figures_json else [],
                )
            )

        return ProvenanceLog(
            session_id=sess_id,
            started_at=started_at,
            input_data=InputData(
                filename=filename or "",
                sha256=sha256 or "",
                n_obs=n_obs or 0,
                n_vars=n_vars or 0,
            ),
            steps=steps,
        )

    async def session_exists(self, session_id: str) -> bool:
        async with self._conn().execute(
            "SELECT 1 FROM sessions WHERE session_id=?", (session_id,)
        ) as cur:
            return await cur.fetchone() is not None

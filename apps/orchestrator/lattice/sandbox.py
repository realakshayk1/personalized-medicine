"""Sandbox abstraction for primitive execution.

Two implementations:
- FakeSandbox: runs primitives in-process via lattice_primitives.run_primitive.
  Safe for local dev because all code goes through the schema-validated,
  version-pinned primitive registry (ADR-001). Used as default for v0.1.
- E2BSandbox: wraps e2b-code-interpreter; used for cloud burst (opt-in).

ADR-003 documents the decision that FakeSandbox is the default for v0.1.
"""

from __future__ import annotations

import io
import json
import os
import tempfile
from typing import Any, Protocol, runtime_checkable

import anndata
from loguru import logger

from lattice.models import PrimitiveResult, SanityWarning

# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------


@runtime_checkable
class SandboxProtocol(Protocol):
    """Abstract sandbox interface."""

    async def upload_anndata(self, adata: anndata.AnnData) -> None:
        """Upload an AnnData object into the sandbox."""
        ...

    async def run_primitive(self, name: str, params: dict[str, Any]) -> PrimitiveResult:
        """Run a named primitive and return the result."""
        ...

    async def download_anndata(self) -> anndata.AnnData:
        """Download the current AnnData from the sandbox."""
        ...

    async def close(self) -> None:
        """Release sandbox resources."""
        ...


# ---------------------------------------------------------------------------
# FakeSandbox — in-process execution (default for v0.1 local-first)
# ---------------------------------------------------------------------------


class FakeSandbox:
    """Executes primitives in-process via lattice_primitives.run_primitive.

    This is architecturally safe for v0.1 because:
    1. All execution flows through the registered, schema-validated primitive
       registry (ADR-001). There is no exec/eval path.
    2. Primitives are version-pinned, human-written, and audited.
    3. The FakeSandbox is never used for arbitrary code execution.

    See ADR-003 in docs/DECISIONS.md.
    """

    def __init__(self) -> None:
        self._adata: anndata.AnnData | None = None
        self._ensure_primitives_registered()

    @staticmethod
    def _ensure_primitives_registered() -> None:
        """Import primitive modules so @primitive decorators fire and register them."""
        import lattice_primitives.all_primitives  # noqa: F401

    async def upload_anndata(self, adata: anndata.AnnData) -> None:
        self._adata = adata
        logger.debug(f"FakeSandbox: uploaded AnnData {adata.n_obs}x{adata.n_vars}")

    async def run_primitive(self, name: str, params: dict[str, Any]) -> PrimitiveResult:
        from lattice_primitives import run_primitive as _run

        if self._adata is None:
            raise RuntimeError("No AnnData loaded — call upload_anndata first")

        logger.info(f"FakeSandbox: running primitive '{name}' params={params}")
        adata_out, warnings, prim_result = _run(name, self._adata, **params)
        self._adata = adata_out

        # Convert lattice_primitives.Warning objects to our Pydantic SanityWarning
        sanity_warnings = [
            SanityWarning(
                severity=w.severity,
                message=w.message,
                primitive=w.primitive,
            )
            for w in warnings
        ]

        return PrimitiveResult(
            primitive=prim_result.primitive,
            params=prim_result.params,
            input_hash=prim_result.input_hash,
            output_hash=prim_result.output_hash,
            duration_sec=prim_result.duration_sec,
            warnings=sanity_warnings,
            user_explanation=prim_result.user_explanation,
            methods_paragraph=prim_result.methods_paragraph,
            n_obs_before=prim_result.n_obs_before,
            n_obs_after=prim_result.n_obs_after,
            n_vars_before=prim_result.n_vars_before,
            n_vars_after=prim_result.n_vars_after,
            figures=list(prim_result.figures),
        )

    async def download_anndata(self) -> anndata.AnnData:
        if self._adata is None:
            raise RuntimeError("No AnnData available")
        return self._adata

    async def close(self) -> None:
        self._adata = None
        logger.debug("FakeSandbox: closed")


# ---------------------------------------------------------------------------
# E2BSandbox — cloud burst (opt-in, not used in v0.1 tests)
# ---------------------------------------------------------------------------


class E2BSandbox:
    """Executes primitives inside an E2B sandbox (cloud burst mode).

    Requires E2B_API_KEY environment variable. Not used by default in v0.1.
    See ADR-003.
    """

    def __init__(self) -> None:
        # Lazy import so missing e2b-code-interpreter doesn't break FakeSandbox users
        try:
            from e2b_code_interpreter import Sandbox

            self._Sandbox = Sandbox
        except ImportError as exc:
            raise ImportError(
                "e2b-code-interpreter is required for E2BSandbox. "
                "Install it or use FakeSandbox for local dev."
            ) from exc
        self._sbx: Any = None
        self._remote_path = "/workspace/data.h5ad"

    async def upload_anndata(self, adata: anndata.AnnData) -> None:
        if self._sbx is None:
            self._sbx = self._Sandbox()
        buf = io.BytesIO()
        adata.write_h5ad(buf)
        buf.seek(0)
        self._sbx.files.write(self._remote_path, buf.read())
        logger.info(f"E2BSandbox: uploaded AnnData {adata.n_obs}x{adata.n_vars}")

    # Sentinels wrap the JSON payload in stdout so we can extract it robustly
    # even if the primitive (or its dependencies) print extra lines.
    _SENTINEL_START = "__LATTICE_RESULT__"
    _SENTINEL_END = "__END__"

    async def run_primitive(self, name: str, params: dict[str, Any]) -> PrimitiveResult:
        if self._sbx is None:
            raise RuntimeError("E2BSandbox not initialized")

        params_json = json.dumps(params)
        code = (
            f"import json\n"
            f"import anndata as ad\n"
            f"from lattice_primitives import run_primitive, PRIMITIVE_REGISTRY\n"
            f"import lattice_primitives.all_primitives\n"
            f"adata = ad.read_h5ad('{self._remote_path}')\n"
            f"params = json.loads('''{params_json}''')\n"
            f"adata_out, warnings, result = run_primitive('{name}', adata, **params)\n"
            f"adata_out.write_h5ad('{self._remote_path}')\n"
            # Wrap the JSON in sentinels so extra stdout lines don't corrupt parsing.
            f"print('{self._SENTINEL_START}' + result.model_dump_json() + '{self._SENTINEL_END}')\n"
        )
        execution = self._sbx.run_code(code)
        return self._parse_execution(execution, name)

    @classmethod
    def _parse_execution(cls, execution: Any, name: str) -> PrimitiveResult:
        """Defensively parse an E2B execution result into a PrimitiveResult.

        Raises RuntimeError on remote error, empty stdout, or malformed JSON.
        Figures (PrimitiveResult.figures) ride along inside the JSON payload.
        """
        # 1. Remote execution error (exception / traceback in the sandbox).
        error = getattr(execution, "error", None)
        if error:
            traceback = (
                getattr(error, "traceback", None)
                or getattr(error, "value", None)
                or str(error)
            )
            raise RuntimeError(f"E2B primitive '{name}' failed: {traceback}")

        # 2. Pull stdout. e2b exposes it as `.text` (logs.stdout is also possible).
        stdout = getattr(execution, "text", None)
        if not stdout:
            logs = getattr(execution, "logs", None)
            stdout_lines = getattr(logs, "stdout", None) if logs is not None else None
            if stdout_lines:
                stdout = "".join(stdout_lines)

        if not stdout or not str(stdout).strip():
            raise RuntimeError(
                f"E2B primitive '{name}' produced no stdout; cannot parse result."
            )

        stdout = str(stdout)

        # 3. Extract the sentinel-wrapped JSON payload.
        start = stdout.find(cls._SENTINEL_START)
        end = stdout.rfind(cls._SENTINEL_END)
        if start != -1 and end != -1 and end > start:
            payload = stdout[start + len(cls._SENTINEL_START) : end]
        else:
            # Fallback: take the last non-empty line and hope it is the JSON.
            non_empty = [ln for ln in stdout.splitlines() if ln.strip()]
            if not non_empty:
                raise RuntimeError(
                    f"E2B primitive '{name}' stdout had no parseable content. "
                    f"Raw stdout: {stdout[:500]!r}"
                )
            payload = non_empty[-1]

        # 4. Parse JSON defensively.
        try:
            result_dict = json.loads(payload)
        except (json.JSONDecodeError, ValueError) as exc:
            raise RuntimeError(
                f"E2B primitive '{name}' returned non-JSON output. "
                f"Raw stdout snippet: {stdout[:500]!r}"
            ) from exc

        # 5. Build the validated Pydantic model (figures included in the dict).
        try:
            return PrimitiveResult(**result_dict)
        except Exception as exc:
            raise RuntimeError(
                f"E2B primitive '{name}' returned JSON that is not a valid "
                f"PrimitiveResult: {exc}. Raw payload snippet: {payload[:500]!r}"
            ) from exc

    async def download_anndata(self) -> anndata.AnnData:
        if self._sbx is None:
            raise RuntimeError("E2BSandbox not initialized")
        data = self._sbx.files.read(self._remote_path)
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
            f.write(data)
            f.flush()
            return anndata.read_h5ad(f.name)

    async def close(self) -> None:
        if self._sbx is not None:
            self._sbx.kill()
            self._sbx = None
        logger.debug("E2BSandbox: closed")


# ---------------------------------------------------------------------------
# Factory — selects sandbox implementation (FakeSandbox is the default; ADR-003)
# ---------------------------------------------------------------------------


def make_sandbox(mode: str | None = None) -> SandboxProtocol:
    """Return a sandbox implementation.

    Selection order:
    1. Explicit ``mode`` argument ("fake" or "e2b").
    2. ``LATTICE_SANDBOX`` environment variable (same values).
    3. Default: FakeSandbox (ADR-003 — in-process is the v0.1 default).

    E2BSandbox is only constructed when explicitly requested. If the optional
    ``e2b-code-interpreter`` dependency is missing, E2BSandbox's constructor
    raises a clear ImportError; that propagates here unchanged.
    """
    resolved = (mode or os.environ.get("LATTICE_SANDBOX") or "fake").strip().lower()

    if resolved == "e2b":
        logger.info("make_sandbox: selecting E2BSandbox (cloud burst)")
        return E2BSandbox()

    if resolved not in ("fake", ""):
        logger.warning(
            f"make_sandbox: unknown sandbox mode {resolved!r}; defaulting to FakeSandbox"
        )

    logger.debug("make_sandbox: selecting FakeSandbox (default, ADR-003)")
    return FakeSandbox()

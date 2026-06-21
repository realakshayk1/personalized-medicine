"""Tests for the sandbox abstraction.

Covers:
- make_sandbox() factory selection (FakeSandbox default; E2B opt-in).
- E2BSandbox.run_primitive output PARSING in isolation (no real E2B needed):
  the sandbox object is faked so we exercise the sentinel-extraction and
  defensive JSON-parsing logic directly.
- FakeSandbox end-to-end (upload AnnData, run a primitive, get PrimitiveResult).

Per ADR-003, FakeSandbox is the default and the only path required for tests.
No E2B credentials are needed.
"""

from __future__ import annotations

import importlib.util
import json
from typing import Any

import anndata
import pytest
from lattice.models import PrimitiveResult
from lattice.sandbox import (
    E2BSandbox,
    FakeSandbox,
    SandboxProtocol,
    make_sandbox,
)

_E2B_INSTALLED = importlib.util.find_spec("e2b_code_interpreter") is not None


# ---------------------------------------------------------------------------
# Fakes for E2B execution-result parsing
# ---------------------------------------------------------------------------


class _FakeError:
    """Mimics the e2b execution `.error` object."""

    def __init__(self, traceback: str) -> None:
        self.traceback = traceback
        self.value = traceback
        self.name = "RuntimeError"

    def __str__(self) -> str:
        return self.value


class _FakeExecution:
    """Mimics the object returned by e2b Sandbox.run_code()."""

    def __init__(self, text: str | None = None, error: Any = None) -> None:
        self.text = text
        self.error = error


def _wrap(payload_json: str) -> str:
    return f"{E2BSandbox._SENTINEL_START}{payload_json}{E2BSandbox._SENTINEL_END}"


def _valid_result_json(with_figure: bool = False) -> str:
    result = PrimitiveResult(
        primitive="calculate_qc_metrics",
        params={"mt_prefix": "MT-"},
        input_hash="abc123",
        output_hash="def456",
        duration_sec=0.05,
        user_explanation="Computed QC metrics.",
        figures=["data:image/png;base64,IVBORw0KGgo="] if with_figure else [],
    )
    return result.model_dump_json()


# ---------------------------------------------------------------------------
# make_sandbox() factory
# ---------------------------------------------------------------------------


class TestMakeSandbox:
    def test_default_is_fake(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("LATTICE_SANDBOX", raising=False)
        sb = make_sandbox()
        assert isinstance(sb, FakeSandbox)
        assert isinstance(sb, SandboxProtocol)

    def test_mode_fake(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("LATTICE_SANDBOX", raising=False)
        assert isinstance(make_sandbox(mode="fake"), FakeSandbox)

    def test_unknown_mode_falls_back_to_fake(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("LATTICE_SANDBOX", raising=False)
        assert isinstance(make_sandbox(mode="bogus"), FakeSandbox)

    def test_env_var_fake(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LATTICE_SANDBOX", "fake")
        assert isinstance(make_sandbox(), FakeSandbox)

    def test_explicit_mode_overrides_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # env says e2b, but explicit mode="fake" must win
        monkeypatch.setenv("LATTICE_SANDBOX", "e2b")
        assert isinstance(make_sandbox(mode="fake"), FakeSandbox)

    @pytest.mark.skipif(
        _E2B_INSTALLED, reason="e2b-code-interpreter is installed; ImportError path n/a"
    )
    def test_mode_e2b_raises_importerror_when_missing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("LATTICE_SANDBOX", raising=False)
        with pytest.raises(ImportError, match="e2b-code-interpreter"):
            make_sandbox(mode="e2b")

    @pytest.mark.skipif(
        not _E2B_INSTALLED, reason="requires e2b-code-interpreter installed"
    )
    def test_mode_e2b_returns_e2b_when_installed(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.delenv("LATTICE_SANDBOX", raising=False)
        sb = make_sandbox(mode="e2b")
        assert isinstance(sb, E2BSandbox)

    @pytest.mark.skipif(
        not _E2B_INSTALLED, reason="requires e2b-code-interpreter installed"
    )
    def test_env_var_e2b_returns_e2b_when_installed(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("LATTICE_SANDBOX", "e2b")
        assert isinstance(make_sandbox(), E2BSandbox)


# ---------------------------------------------------------------------------
# E2BSandbox output parsing (in isolation — no real E2B)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not _E2B_INSTALLED, reason="requires e2b-code-interpreter installed"
)
class TestE2BParsing:
    def _make_sandbox(self) -> E2BSandbox:
        return E2BSandbox()

    def test_parses_sentinel_wrapped_json(self) -> None:
        execution = _FakeExecution(text=_wrap(_valid_result_json()))
        result = E2BSandbox._parse_execution(execution, "calculate_qc_metrics")
        assert isinstance(result, PrimitiveResult)
        assert result.primitive == "calculate_qc_metrics"
        assert result.output_hash == "def456"

    def test_figures_survive_roundtrip(self) -> None:
        execution = _FakeExecution(text=_wrap(_valid_result_json(with_figure=True)))
        result = E2BSandbox._parse_execution(execution, "plot_something")
        assert len(result.figures) == 1
        assert result.figures[0].startswith("data:image/png;base64,")

    def test_extra_stdout_lines_are_ignored(self) -> None:
        noisy = (
            "WARNING: some library chatter\n"
            "progress: 50%\n"
            + _wrap(_valid_result_json())
            + "\ntrailing noise after sentinel\n"
        )
        execution = _FakeExecution(text=noisy)
        result = E2BSandbox._parse_execution(execution, "calculate_qc_metrics")
        assert result.primitive == "calculate_qc_metrics"

    def test_fallback_last_line_when_no_sentinel(self) -> None:
        # No sentinels; last non-empty line is valid JSON.
        text = "some log line\n" + _valid_result_json() + "\n"
        execution = _FakeExecution(text=text)
        result = E2BSandbox._parse_execution(execution, "calculate_qc_metrics")
        assert result.primitive == "calculate_qc_metrics"

    def test_execution_error_raises_runtimeerror(self) -> None:
        execution = _FakeExecution(
            text=None,
            error=_FakeError("Traceback (most recent call last): ValueError: boom"),
        )
        with pytest.raises(RuntimeError, match="boom"):
            E2BSandbox._parse_execution(execution, "calculate_qc_metrics")

    def test_empty_stdout_raises_runtimeerror(self) -> None:
        execution = _FakeExecution(text="")
        with pytest.raises(RuntimeError, match="no stdout"):
            E2BSandbox._parse_execution(execution, "calculate_qc_metrics")

    def test_whitespace_only_stdout_raises(self) -> None:
        execution = _FakeExecution(text="   \n  \n")
        with pytest.raises(RuntimeError, match="no stdout"):
            E2BSandbox._parse_execution(execution, "calculate_qc_metrics")

    def test_malformed_json_raises_with_snippet(self) -> None:
        execution = _FakeExecution(text=_wrap("{not valid json at all"))
        with pytest.raises(RuntimeError, match="non-JSON output"):
            E2BSandbox._parse_execution(execution, "calculate_qc_metrics")

    def test_valid_json_wrong_shape_raises(self) -> None:
        # Valid JSON, but not a PrimitiveResult (missing required fields).
        execution = _FakeExecution(text=_wrap(json.dumps({"foo": "bar"})))
        with pytest.raises(RuntimeError, match="not a valid"):
            E2BSandbox._parse_execution(execution, "calculate_qc_metrics")

    def test_run_primitive_uses_parser(self) -> None:
        """End-to-end of run_primitive with a faked _sbx whose run_code
        returns sentinel-wrapped JSON. No real E2B involved."""
        import asyncio

        class _FakeSbx:
            def __init__(self) -> None:
                self.last_code: str | None = None

            def run_code(self, code: str) -> _FakeExecution:
                self.last_code = code
                return _FakeExecution(text=_wrap(_valid_result_json()))

        sb = self._make_sandbox()
        sb._sbx = _FakeSbx()  # type: ignore[assignment]
        result = asyncio.run(
            sb.run_primitive("calculate_qc_metrics", {"mt_prefix": "MT-"})
        )
        assert isinstance(result, PrimitiveResult)
        assert result.primitive == "calculate_qc_metrics"
        # Remote code must still import the full primitive registry.
        assert "import lattice_primitives.all_primitives" in sb._sbx.last_code  # type: ignore[union-attr]


# ---------------------------------------------------------------------------
# FakeSandbox end-to-end
# ---------------------------------------------------------------------------


class TestFakeSandboxEndToEnd:
    @pytest.mark.asyncio
    async def test_upload_run_download(self, synthetic_adata: anndata.AnnData) -> None:
        sb = FakeSandbox()
        await sb.upload_anndata(synthetic_adata)
        result = await sb.run_primitive("calculate_qc_metrics", {"mt_prefix": "MT-"})
        assert isinstance(result, PrimitiveResult)
        assert result.primitive == "calculate_qc_metrics"
        assert result.input_hash
        assert result.output_hash
        adata_out = await sb.download_anndata()
        assert adata_out.n_obs == synthetic_adata.n_obs
        await sb.close()

    @pytest.mark.asyncio
    async def test_run_without_upload_raises(self) -> None:
        sb = FakeSandbox()
        with pytest.raises(RuntimeError, match="No AnnData"):
            await sb.run_primitive("calculate_qc_metrics", {"mt_prefix": "MT-"})

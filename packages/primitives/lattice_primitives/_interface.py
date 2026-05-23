"""Frozen primitive interface contract.

THIS FILE IS A FROZEN CONTRACT. Subagents in Wave 2 MUST NOT modify it.
If you believe the contract is wrong, emit `BLOCKED: contract change requested
in _interface.py` and stop. Do not silently edit.

Imports allowed from this module:
    Primitive, PrimitiveResult, ParamSpec, AnnDataState, Severity, Warning,
    PRIMITIVE_REGISTRY, primitive, sanity_check

Implementation of the @primitive and @sanity_check decorators, the registry
dict, and schema validation lives in `registry.py` (Wave 2-A). This module
defines only the Pydantic models and protocol that all consumers can rely on.
"""

from __future__ import annotations

from typing import Any, Callable, Literal, Protocol

from pydantic import BaseModel, ConfigDict, Field

Severity = Literal["info", "warn", "error"]


class Warning_(BaseModel):
    """A sanity-check finding emitted by a primitive after execution."""

    model_config = ConfigDict(frozen=True)
    severity: Severity
    message: str
    primitive: str


# Public alias matching AGENTS.md prose; underscore avoids shadowing builtin.
Warning = Warning_


class ParamSpec(BaseModel):
    """Schema for one parameter on a primitive."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    type_: str = Field(alias="type", description="Python type name as string: int, float, bool, str, list[str], etc.")
    default: Any = None
    description: str = ""
    choices: list[Any] | None = None
    min_value: float | None = None
    max_value: float | None = None


class AnnDataState(BaseModel):
    """Schema describing required or produced AnnData state.

    `obs_columns` / `var_columns`: column names that must be present (requires)
    or are added (produces). Required columns are validated by the registry
    BEFORE primitive execution; produced columns are validated AFTER.

    `layer_state` is a string tag tracking the data's transformation lineage
    (e.g. "raw_counts", "log_normalized"). The registry refuses to run a
    primitive whose `requires.layer_state` does not match the current state.

    `adds_layers` / `adds_obsm` / `adds_uns` declare additions to those slots.
    """

    model_config = ConfigDict(frozen=True)
    obs_columns: list[str] = Field(default_factory=list)
    var_columns: list[str] = Field(default_factory=list)
    layer_state: str | None = None
    adds_layers: list[str] = Field(default_factory=list)
    adds_obsm: list[str] = Field(default_factory=list)
    adds_uns: list[str] = Field(default_factory=list)
    # Columns whose VALUES must be preserved bit-for-bit through this primitive.
    # The canonical use is `["condition"]` — protects against case/control swap.
    immutable_obs_columns: list[str] = Field(default_factory=list)


class PrimitiveResult(BaseModel):
    """Structured result of one primitive invocation.

    The AnnData itself is not embedded here — it lives in the sandbox. This
    is the JSON-serializable record streamed to the frontend and persisted
    in the provenance log.
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    primitive: str
    params: dict[str, Any]
    input_hash: str
    output_hash: str
    duration_sec: float
    warnings: list[Warning_] = Field(default_factory=list)
    user_explanation: str = ""
    methods_paragraph: str = ""
    n_obs_before: int | None = None
    n_obs_after: int | None = None
    n_vars_before: int | None = None
    n_vars_after: int | None = None
    figures: list[str] = Field(default_factory=list)  # paths or data-URIs


class Primitive(BaseModel):
    """Registry entry for one primitive. Populated by the @primitive decorator."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    name: str
    category: Literal["qc", "preprocess", "dim_reduce", "cluster", "integrate", "annotate", "de", "plot"]
    requires: AnnDataState
    produces: AnnDataState
    params: dict[str, ParamSpec]
    scanpy_version: str
    citation_key: str
    func: Callable[..., Any]
    sanity_check_func: Callable[..., list[Warning_]] | None = None
    template_path: str | None = None


# The registry dict is populated by @primitive at import time.
# Subagent A owns the implementation in registry.py and re-exports this name.
PRIMITIVE_REGISTRY: dict[str, Primitive] = {}


class PrimitiveCallable(Protocol):
    """Callable shape for a registered primitive (post-decoration)."""

    __primitive_name__: str

    def __call__(self, adata: Any, **params: Any) -> Any: ...


# Decorator stubs — real implementations live in registry.py. These exist
# only so type-checkers can resolve imports from this module during Wave 1.
# Subagent A MUST replace these with the real implementations in registry.py
# and re-export from this module's __init__ alongside (NOT inside this file).
def primitive(**kwargs: Any) -> Callable[[Callable[..., Any]], Callable[..., Any]]:  # pragma: no cover
    raise NotImplementedError("Implemented in registry.py — Wave 2-A")


def sanity_check(primitive_name: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:  # pragma: no cover
    raise NotImplementedError("Implemented in registry.py — Wave 2-A")

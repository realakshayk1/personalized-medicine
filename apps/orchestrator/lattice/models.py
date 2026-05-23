"""Pydantic models mirroring packages/sdk/src/api.ts exactly.

Field names MUST match the TypeScript interface field names one-to-one.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Primitive + Plan shapes
# ---------------------------------------------------------------------------


class SanityWarning(BaseModel):
    severity: Literal["info", "warn", "error"]
    message: str
    primitive: str


class PrimitiveInvocation(BaseModel):
    primitive: str
    params: dict[str, Any] = Field(default_factory=dict)
    rationale: str | None = None


class Plan(BaseModel):
    session_id: str
    steps: list[PrimitiveInvocation]
    rationale: str
    workflow: str


class PrimitiveResult(BaseModel):
    primitive: str
    params: dict[str, Any]
    input_hash: str
    output_hash: str
    duration_sec: float
    warnings: list[SanityWarning] = Field(default_factory=list)
    user_explanation: str = ""
    methods_paragraph: str = ""
    n_obs_before: int | None = None
    n_obs_after: int | None = None
    n_vars_before: int | None = None
    n_vars_after: int | None = None
    figures: list[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Session shapes
# ---------------------------------------------------------------------------


class SessionCreateResponse(BaseModel):
    session_id: str


class UploadResponse(BaseModel):
    session_id: str
    filename: str
    sha256: str
    n_obs: int
    n_vars: int
    obs_columns: list[str]


class PlanRequest(BaseModel):
    user_message: str


class ExecuteRequest(BaseModel):
    plan: Plan


# ---------------------------------------------------------------------------
# Streaming event shapes
# ---------------------------------------------------------------------------


class StepStartedEvent(BaseModel):
    type: Literal["step_started"] = "step_started"
    step_index: int
    primitive: str


class StepCompletedEvent(BaseModel):
    type: Literal["step_completed"] = "step_completed"
    step_index: int
    result: PrimitiveResult


class StepFailedEvent(BaseModel):
    type: Literal["step_failed"] = "step_failed"
    step_index: int
    primitive: str
    error: str


class PlanCompletedEvent(BaseModel):
    type: Literal["plan_completed"] = "plan_completed"
    n_steps: int


ExecuteEvent = StepStartedEvent | StepCompletedEvent | StepFailedEvent | PlanCompletedEvent


# ---------------------------------------------------------------------------
# Provenance shapes
# ---------------------------------------------------------------------------


class ProvenanceEntry(PrimitiveResult):
    step_id: int
    started_at: str


class InputData(BaseModel):
    filename: str
    sha256: str
    n_obs: int
    n_vars: int


class ProvenanceLog(BaseModel):
    session_id: str
    started_at: str
    input_data: InputData
    steps: list[ProvenanceEntry]

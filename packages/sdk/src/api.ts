/**
 * FROZEN API CONTRACT. Subagents in Wave 2 MUST NOT modify this file.
 * If the contract is wrong, emit `BLOCKED: contract change requested in api.ts`
 * and stop. Do not silently edit.
 *
 * This file is the single source of truth for the wire shape between the
 * desktop frontend and the orchestrator. Both sides import from here.
 */

export type Severity = "info" | "warn" | "error";

export interface SanityWarning {
  severity: Severity;
  message: string;
  primitive: string;
}

export interface PrimitiveInvocation {
  primitive: string;
  params: Record<string, unknown>;
  rationale?: string;
}

export interface Plan {
  session_id: string;
  steps: PrimitiveInvocation[];
  rationale: string;
  workflow: string;
}

export interface PrimitiveResult {
  primitive: string;
  params: Record<string, unknown>;
  input_hash: string;
  output_hash: string;
  duration_sec: number;
  warnings: SanityWarning[];
  user_explanation: string;
  methods_paragraph: string;
  n_obs_before: number | null;
  n_obs_after: number | null;
  n_vars_before: number | null;
  n_vars_after: number | null;
  figures: string[];
}

export interface SessionCreateResponse {
  session_id: string;
}

export interface UploadResponse {
  session_id: string;
  filename: string;
  sha256: string;
  n_obs: number;
  n_vars: number;
  obs_columns: string[];
}

export interface PlanRequest {
  user_message: string;
}

export interface ExecuteRequest {
  plan: Plan;
}

// Server-sent events from /sessions/{id}/execute
export type ExecuteEvent =
  | { type: "step_started"; step_index: number; primitive: string }
  | { type: "step_completed"; step_index: number; result: PrimitiveResult }
  | { type: "step_failed"; step_index: number; primitive: string; error: string }
  | { type: "plan_completed"; n_steps: number };

export interface ProvenanceEntry extends PrimitiveResult {
  step_id: number;
  started_at: string;
}

export interface ProvenanceLog {
  session_id: string;
  started_at: string;
  input_data: {
    filename: string;
    sha256: string;
    n_obs: number;
    n_vars: number;
  };
  steps: ProvenanceEntry[];
}

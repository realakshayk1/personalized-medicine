"use client";

import { useState } from "react";
import { CheckCircle2, Circle, XCircle, Loader2, Play } from "lucide-react";
import { cn } from "@/lib/utils";
import { executePlan } from "@/lib/api";
import type { Plan, PrimitiveResult, ExecuteEvent } from "@lattice/sdk/api";

interface PlanPanelProps {
  plan: Plan | null;
  sessionId: string;
  onStepResult: (result: PrimitiveResult) => void;
}

type StepStatus = "pending" | "running" | "done" | "error";

interface StepState {
  status: StepStatus;
  result?: PrimitiveResult;
  error?: string;
}

export function PlanPanel({ plan, sessionId, onStepResult }: PlanPanelProps) {
  const [stepStates, setStepStates] = useState<StepState[]>([]);
  const [executing, setExecuting] = useState(false);
  const [executionDone, setExecutionDone] = useState(false);

  const runPlan = async () => {
    if (!plan || executing) return;
    setExecuting(true);
    setExecutionDone(false);
    setStepStates(plan.steps.map(() => ({ status: "pending" })));

    try {
      await executePlan(sessionId, plan, (event: ExecuteEvent) => {
        if (event.type === "step_started") {
          setStepStates((prev) => {
            const next = [...prev];
            if (next[event.step_index]) {
              next[event.step_index] = { status: "running" };
            }
            return next;
          });
        } else if (event.type === "step_completed") {
          onStepResult(event.result);
          setStepStates((prev) => {
            const next = [...prev];
            if (next[event.step_index]) {
              next[event.step_index] = { status: "done", result: event.result };
            }
            return next;
          });
        } else if (event.type === "step_failed") {
          setStepStates((prev) => {
            const next = [...prev];
            if (next[event.step_index]) {
              next[event.step_index] = {
                status: "error",
                error: event.error,
              };
            }
            return next;
          });
        } else if (event.type === "plan_completed") {
          setExecutionDone(true);
        }
      });
    } catch (err) {
      console.error("executePlan error:", err);
    } finally {
      setExecuting(false);
    }
  };

  if (!plan) {
    return (
      <div className="flex items-center justify-center h-full text-xs text-[hsl(215,20%,45%)] px-4 text-center">
        No plan yet. Describe your analysis in the chat to generate one.
      </div>
    );
  }

  const stepsWithState = plan.steps.map((step, i) => ({
    step,
    state: stepStates[i] ?? { status: "pending" as StepStatus },
  }));

  return (
    <div className="flex flex-col h-full p-3 gap-3 overflow-hidden">
      <div className="flex items-center justify-between shrink-0">
        <div className="text-[11px] text-[hsl(215,20%,55%)]">
          {plan.workflow} — {plan.steps.length} steps
        </div>
        <button
          onClick={() => void runPlan()}
          disabled={executing || executionDone}
          className={cn(
            "flex items-center gap-1.5 px-3 py-1 rounded text-[11px] font-medium transition-colors",
            "bg-[hsl(217,91%,60%)] text-white",
            "disabled:opacity-40 hover:enabled:bg-[hsl(217,91%,55%)]"
          )}
          data-testid="run-plan-button"
        >
          {executing ? (
            <Loader2 className="w-3 h-3 animate-spin" />
          ) : (
            <Play className="w-3 h-3" />
          )}
          {executing ? "Running…" : executionDone ? "Done" : "Run Plan"}
        </button>
      </div>

      {plan.rationale && (
        <p className="text-[11px] text-[hsl(215,20%,55%)] leading-relaxed shrink-0">
          {plan.rationale}
        </p>
      )}

      <div className="flex-1 overflow-y-auto space-y-2">
        {stepsWithState.map(({ step, state }, i) => (
          <div
            key={i}
            className={cn(
              "rounded border p-3 transition-colors",
              state.status === "running" &&
                "border-[hsl(217,91%,60%,0.5)] bg-[hsl(217,91%,60%,0.05)]",
              state.status === "done" && "border-[hsl(217,32%,18%)]",
              state.status === "error" && "border-red-500/30",
              state.status === "pending" && "border-[hsl(217,32%,14%)]"
            )}
          >
            <div className="flex items-start gap-2">
              <span className="mt-0.5 shrink-0">
                {state.status === "pending" && (
                  <Circle className="w-3.5 h-3.5 text-[hsl(215,20%,35%)]" />
                )}
                {state.status === "running" && (
                  <Loader2 className="w-3.5 h-3.5 text-[hsl(217,91%,60%)] animate-spin" />
                )}
                {state.status === "done" && (
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400" />
                )}
                {state.status === "error" && (
                  <XCircle className="w-3.5 h-3.5 text-red-400" />
                )}
              </span>
              <div className="flex-1 min-w-0">
                <p className="text-xs font-mono text-[hsl(210,40%,88%)]">
                  {step.primitive}
                </p>
                <div className="mt-1 flex flex-wrap gap-x-3 gap-y-0.5">
                  {Object.entries(step.params).map(([k, v]) => (
                    <span
                      key={k}
                      className="text-[10px] font-mono text-[hsl(215,20%,55%)]"
                    >
                      {k}={JSON.stringify(v)}
                    </span>
                  ))}
                </div>
                {step.rationale && (
                  <p className="mt-1 text-[10px] text-[hsl(215,20%,45%)]">
                    {step.rationale}
                  </p>
                )}
                {state.result?.user_explanation && (
                  <p className="mt-1.5 text-[11px] text-emerald-300/80">
                    {state.result.user_explanation}
                  </p>
                )}
                {state.error && (
                  <p className="mt-1 text-[11px] text-red-400">{state.error}</p>
                )}
                {state.result?.warnings && state.result.warnings.length > 0 && (
                  <div className="mt-1 space-y-0.5">
                    {state.result.warnings.map((w, wi) => (
                      <p
                        key={wi}
                        className={cn(
                          "text-[10px]",
                          w.severity === "error"
                            ? "text-red-400"
                            : w.severity === "warn"
                            ? "text-yellow-400"
                            : "text-[hsl(215,20%,55%)]"
                        )}
                      >
                        [{w.severity}] {w.message}
                      </p>
                    ))}
                  </div>
                )}
              </div>
              <span className="text-[10px] text-[hsl(215,20%,40%)] shrink-0 font-mono">
                {(i + 1).toString().padStart(2, "0")}
              </span>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

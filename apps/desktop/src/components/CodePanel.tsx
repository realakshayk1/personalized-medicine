"use client";

import { Code2 } from "lucide-react";
import type { Plan } from "@lattice/sdk/api";

interface CodePanelProps {
  plan: Plan | null;
}

function planToPseudoPython(plan: Plan): string {
  const lines: string[] = [
    "import scanpy as sc",
    "import anndata as ad",
    "",
    `# Session: ${plan.session_id}`,
    `# Workflow: ${plan.workflow}`,
    "",
    "adata = ad.read_h5ad('data.h5ad')",
    "",
  ];

  for (const step of plan.steps) {
    lines.push(`# ${step.rationale ?? step.primitive}`);
    const paramStr = Object.entries(step.params)
      .map(([k, v]) => `${k}=${JSON.stringify(v)}`)
      .join(", ");
    lines.push(`primitives.${step.primitive}(adata, ${paramStr})`);
    lines.push("");
  }

  return lines.join("\n");
}

export function CodePanel({ plan }: CodePanelProps) {
  return (
    <div className="flex flex-col h-full" data-testid="code-panel">
      {plan === null ? (
        <div className="flex flex-col items-center justify-center flex-1 gap-3 text-center px-6">
          <Code2 className="w-8 h-8 text-[hsl(217,32%,22%)]" />
          <p className="text-xs text-[hsl(215,20%,40%)]">
            Code preview will appear here once a plan is generated.
          </p>
        </div>
      ) : (
        <div className="flex-1 overflow-auto p-3">
          <pre className="text-[11px] font-mono leading-relaxed text-[hsl(215,20%,72%)] whitespace-pre">
            <code>{planToPseudoPython(plan)}</code>
          </pre>
        </div>
      )}
    </div>
  );
}

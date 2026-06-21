"use client";

import { useState } from "react";
import { Code2, Download, Loader2 } from "lucide-react";
import type { Plan } from "@lattice/sdk/api";
import { exportNotebook, type ExportFormat } from "@/lib/api";

interface CodePanelProps {
  plan: Plan | null;
  sessionId: string;
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

function triggerDownload(filename: string, content: string): void {
  const blob = new Blob([content], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

export function CodePanel({ plan, sessionId }: CodePanelProps) {
  const [exporting, setExporting] = useState<ExportFormat | null>(null);
  const [exportError, setExportError] = useState<string | null>(null);
  const [exported, setExported] = useState<{
    filename: string;
    content: string;
  } | null>(null);

  const runExport = async (format: ExportFormat) => {
    if (exporting) return;
    setExporting(format);
    setExportError(null);
    try {
      const result = await exportNotebook(sessionId, format);
      setExported(result);
      triggerDownload(result.filename, result.content);
    } catch (err) {
      setExportError(err instanceof Error ? err.message : String(err));
    } finally {
      setExporting(null);
    }
  };

  return (
    <div className="flex flex-col h-full" data-testid="code-panel">
      <div className="flex items-center justify-between gap-2 px-3 py-2 border-b border-[hsl(217,32%,16%)] shrink-0">
        <span className="text-[11px] text-[hsl(215,20%,55%)]">
          {exported
            ? exported.filename
            : "Reproducible code preview"}
        </span>
        <div className="flex items-center gap-1.5">
          <button
            onClick={() => void runExport("ipynb")}
            disabled={exporting !== null}
            className="flex items-center gap-1.5 px-2.5 py-1 rounded text-[11px] font-medium transition-colors bg-[hsl(217,91%,60%)] text-white disabled:opacity-40 hover:enabled:bg-[hsl(217,91%,55%)]"
            data-testid="export-notebook-button"
          >
            {exporting === "ipynb" ? (
              <Loader2 className="w-3 h-3 animate-spin" />
            ) : (
              <Download className="w-3 h-3" />
            )}
            Export notebook
          </button>
          <button
            onClick={() => void runExport("py")}
            disabled={exporting !== null}
            className="flex items-center gap-1.5 px-2.5 py-1 rounded text-[11px] font-medium transition-colors border border-[hsl(217,32%,22%)] text-[hsl(215,20%,72%)] disabled:opacity-40 hover:enabled:border-[hsl(217,91%,60%)]"
            data-testid="export-py-button"
          >
            {exporting === "py" ? (
              <Loader2 className="w-3 h-3 animate-spin" />
            ) : (
              <Download className="w-3 h-3" />
            )}
            .py
          </button>
        </div>
      </div>

      {exportError && (
        <p className="px-3 py-1.5 text-[11px] text-red-400 shrink-0">
          {exportError}
        </p>
      )}

      {exported ? (
        <div className="flex-1 overflow-auto p-3">
          <pre className="text-[11px] font-mono leading-relaxed text-[hsl(215,20%,72%)] whitespace-pre">
            <code>{exported.content}</code>
          </pre>
        </div>
      ) : plan === null ? (
        <div className="flex flex-col items-center justify-center flex-1 gap-3 text-center px-6">
          <Code2 className="w-8 h-8 text-[hsl(217,32%,22%)]" />
          <p className="text-xs text-[hsl(215,20%,40%)]">
            Code preview will appear here once a plan is generated. Run the plan,
            then export a runnable notebook.
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

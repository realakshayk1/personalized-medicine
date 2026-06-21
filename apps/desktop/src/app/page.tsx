"use client";

import { useState } from "react";
import {
  PanelGroup,
  Panel,
  PanelResizeHandle,
} from "react-resizable-panels";
import { ChatPanel } from "@/components/ChatPanel";
import { FigurePanel } from "@/components/FigurePanel";
import { PlanPanel } from "@/components/PlanPanel";
import { CodePanel } from "@/components/CodePanel";
import { FileDrop } from "@/components/FileDrop";
import type { Plan, UploadResponse, PrimitiveResult } from "@lattice/sdk/api";

type ActiveTab = "plan" | "code";

export default function Home() {
  const [uploadResponse, setUploadResponse] = useState<UploadResponse | null>(
    null
  );
  const [plan, setPlan] = useState<Plan | null>(null);
  const [latestResult, setLatestResult] = useState<PrimitiveResult | null>(
    null
  );
  const [activeTab, setActiveTab] = useState<ActiveTab>("plan");

  const sessionId = uploadResponse?.session_id ?? "demo-session";

  return (
    <div className="h-full flex flex-col bg-[hsl(222,47%,6%)]">
      {/* Top bar */}
      <header className="flex items-center px-4 py-2 border-b border-[hsl(217,32%,16%)] shrink-0">
        <span className="text-sm font-semibold tracking-wide text-[hsl(210,40%,92%)]">
          Lattice
        </span>
        {uploadResponse && (
          <span className="ml-4 text-xs text-[hsl(215,20%,55%)]">
            {uploadResponse.filename} — {uploadResponse.n_obs.toLocaleString()}{" "}
            cells × {uploadResponse.n_vars.toLocaleString()} genes
          </span>
        )}
      </header>

      <div className="flex-1 overflow-hidden">
        <PanelGroup direction="horizontal" className="h-full">
          {/* Left: Chat panel */}
          <Panel
            defaultSize={30}
            minSize={20}
            data-testid="panel-chat"
            className="flex flex-col"
          >
            <div className="flex flex-col h-full p-3 gap-3">
              {!uploadResponse && (
                <FileDrop
                  sessionId={sessionId}
                  onUpload={setUploadResponse}
                />
              )}
              {uploadResponse && (
                <div className="rounded border border-[hsl(217,32%,16%)] bg-[hsl(222,47%,8%)] p-3 text-xs space-y-1">
                  <div className="text-[hsl(215,20%,55%)] font-mono uppercase tracking-wider text-[10px]">
                    Dataset
                  </div>
                  <div className="text-[hsl(210,40%,92%)] font-medium">
                    {uploadResponse.filename}
                  </div>
                  <div className="text-[hsl(215,20%,55%)]">
                    {uploadResponse.n_obs.toLocaleString()} obs ×{" "}
                    {uploadResponse.n_vars.toLocaleString()} vars
                  </div>
                  {uploadResponse.obs_columns.length > 0 && (
                    <div className="text-[hsl(215,20%,55%)]">
                      Columns:{" "}
                      <span className="font-mono">
                        {uploadResponse.obs_columns.slice(0, 6).join(", ")}
                        {uploadResponse.obs_columns.length > 6 && " …"}
                      </span>
                    </div>
                  )}
                </div>
              )}
              <div className="flex-1 overflow-hidden">
                <ChatPanel
                  sessionId={sessionId}
                  onPlanUpdate={setPlan}
                />
              </div>
            </div>
          </Panel>

          <PanelResizeHandle className="w-px bg-[hsl(217,32%,16%)] hover:bg-[hsl(217,91%,60%)] transition-colors" />

          {/* Right side */}
          <Panel defaultSize={70} minSize={40} className="flex flex-col">
            <PanelGroup direction="vertical" className="h-full">
              {/* Top-right: Figure panel */}
              <Panel
                defaultSize={55}
                minSize={20}
                data-testid="panel-figure"
                className="flex flex-col"
              >
                <FigurePanel result={latestResult} />
              </Panel>

              <PanelResizeHandle className="h-px bg-[hsl(217,32%,16%)] hover:bg-[hsl(217,91%,60%)] transition-colors" />

              {/* Bottom-right: Code / Plan tabs */}
              <Panel
                defaultSize={45}
                minSize={20}
                data-testid="panel-bottom-right"
                className="flex flex-col"
              >
                <div className="flex border-b border-[hsl(217,32%,16%)] shrink-0">
                  <button
                    onClick={() => setActiveTab("plan")}
                    className={`px-4 py-2 text-xs font-medium transition-colors ${
                      activeTab === "plan"
                        ? "text-[hsl(210,40%,92%)] border-b-2 border-[hsl(217,91%,60%)]"
                        : "text-[hsl(215,20%,55%)] hover:text-[hsl(210,40%,92%)]"
                    }`}
                  >
                    Plan
                  </button>
                  <button
                    onClick={() => setActiveTab("code")}
                    className={`px-4 py-2 text-xs font-medium transition-colors ${
                      activeTab === "code"
                        ? "text-[hsl(210,40%,92%)] border-b-2 border-[hsl(217,91%,60%)]"
                        : "text-[hsl(215,20%,55%)] hover:text-[hsl(210,40%,92%)]"
                    }`}
                  >
                    Code
                  </button>
                </div>
                <div className="flex-1 overflow-hidden">
                  {activeTab === "plan" ? (
                    <PlanPanel
                      plan={plan}
                      sessionId={sessionId}
                      onStepResult={setLatestResult}
                    />
                  ) : (
                    <CodePanel plan={plan} sessionId={sessionId} />
                  )}
                </div>
              </Panel>
            </PanelGroup>
          </Panel>
        </PanelGroup>
      </div>
    </div>
  );
}

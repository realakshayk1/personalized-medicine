"use client";

import { ImageOff } from "lucide-react";
import type { PrimitiveResult } from "@lattice/sdk/api";

interface FigurePanelProps {
  result: PrimitiveResult | null;
}

export function FigurePanel({ result }: FigurePanelProps) {
  const figures = result?.figures ?? [];

  return (
    <div className="flex flex-col h-full" data-testid="figure-panel">
      <div className="flex items-center px-3 py-2 border-b border-[hsl(217,32%,16%)] shrink-0">
        <span className="text-[11px] font-medium text-[hsl(215,20%,55%)] uppercase tracking-wider">
          Figures
        </span>
        {result && (
          <span className="ml-2 text-[10px] text-[hsl(215,20%,40%)] font-mono">
            {result.primitive}
          </span>
        )}
      </div>

      {figures.length === 0 ? (
        <div className="flex flex-col items-center justify-center flex-1 gap-3 text-center px-6">
          <ImageOff className="w-8 h-8 text-[hsl(217,32%,22%)]" />
          <p className="text-xs text-[hsl(215,20%,40%)]">
            Figures will appear here after execution.
          </p>
          <p className="text-[11px] text-[hsl(215,20%,33%)]">
            Generate a plan and click Run Plan to start.
          </p>
        </div>
      ) : (
        <div className="flex-1 overflow-auto p-3 grid gap-3 auto-rows-min grid-cols-1">
          {figures.map((src, i) => (
            <div
              key={i}
              className="rounded border border-[hsl(217,32%,16%)] overflow-hidden bg-[hsl(222,47%,5%)]"
            >
              {/* reason: figure paths from backend can be URLs or data URIs */}
              {/* eslint-disable-next-line @next/next/no-img-element */}
              <img
                src={src}
                alt={`Figure ${i + 1} from ${result?.primitive ?? "analysis"}`}
                className="w-full h-auto block"
                loading="lazy"
              />
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

"use client";

import { useState, useCallback, useRef } from "react";
import { UploadCloud } from "lucide-react";
import { cn } from "@/lib/utils";
import { uploadH5ad } from "@/lib/api";
import type { UploadResponse } from "@lattice/sdk/api";

interface FileDropProps {
  sessionId: string;
  onUpload: (response: UploadResponse) => void;
}

type DropState = "idle" | "dragging" | "uploading" | "error";

export function FileDrop({ sessionId, onUpload }: FileDropProps) {
  const [dropState, setDropState] = useState<DropState>("idle");
  const [progress, setProgress] = useState(0);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const handleFile = useCallback(
    async (file: File) => {
      if (!file.name.endsWith(".h5ad")) {
        setErrorMsg("Only .h5ad files are supported.");
        setDropState("error");
        return;
      }

      setDropState("uploading");
      setProgress(0);
      setErrorMsg(null);

      try {
        const response = await uploadH5ad(sessionId, file, setProgress);
        onUpload(response);
      } catch (err) {
        setErrorMsg(err instanceof Error ? err.message : "Upload failed.");
        setDropState("error");
      }
    },
    [sessionId, onUpload]
  );

  const onDrop = useCallback(
    (e: React.DragEvent<HTMLDivElement>) => {
      e.preventDefault();
      setDropState("idle");
      const file = e.dataTransfer.files[0];
      if (file) handleFile(file);
    },
    [handleFile]
  );

  const onDragOver = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    setDropState("dragging");
  };

  const onDragLeave = () => {
    setDropState("idle");
  };

  const onInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) handleFile(file);
  };

  return (
    <div
      data-testid="file-drop-zone"
      onDrop={onDrop}
      onDragOver={onDragOver}
      onDragLeave={onDragLeave}
      onClick={() => dropState === "idle" && inputRef.current?.click()}
      className={cn(
        "flex flex-col items-center justify-center gap-2 rounded border border-dashed p-6 text-center cursor-pointer select-none transition-colors",
        dropState === "dragging" &&
          "border-[hsl(217,91%,60%)] bg-[hsl(217,91%,60%,0.08)]",
        dropState === "idle" &&
          "border-[hsl(217,32%,22%)] hover:border-[hsl(217,32%,35%)]",
        dropState === "uploading" && "border-[hsl(217,32%,22%)] cursor-default",
        dropState === "error" && "border-red-500/50"
      )}
      role="button"
      aria-label="Drop .h5ad file here or click to browse"
    >
      <input
        ref={inputRef}
        type="file"
        accept=".h5ad"
        className="hidden"
        onChange={onInputChange}
        aria-hidden="true"
      />

      <UploadCloud
        className={cn(
          "w-7 h-7",
          dropState === "dragging"
            ? "text-[hsl(217,91%,60%)]"
            : "text-[hsl(215,20%,45%)]"
        )}
      />

      {dropState === "uploading" ? (
        <>
          <p className="text-xs text-[hsl(215,20%,55%)]">
            Uploading… {progress}%
          </p>
          <div className="w-full rounded-full h-1 bg-[hsl(217,32%,16%)]">
            <div
              className="h-1 rounded-full bg-[hsl(217,91%,60%)] transition-all"
              style={{ width: `${progress}%` }}
            />
          </div>
        </>
      ) : dropState === "error" ? (
        <p className="text-xs text-red-400">{errorMsg}</p>
      ) : (
        <>
          <p className="text-xs font-medium text-[hsl(210,40%,85%)]">
            Drop <code className="font-mono">.h5ad</code> file here
          </p>
          <p className="text-[11px] text-[hsl(215,20%,50%)]">
            or click to browse
          </p>
        </>
      )}
    </div>
  );
}

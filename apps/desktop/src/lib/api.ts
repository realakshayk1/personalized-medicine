import type {
  SessionCreateResponse,
  UploadResponse,
  Plan,
  PlanRequest,
  ExecuteEvent,
} from "@lattice/sdk/api";

const BASE_URL =
  process.env["NEXT_PUBLIC_ORCHESTRATOR_URL"] ?? "http://localhost:8000";

export async function createSession(): Promise<SessionCreateResponse> {
  const res = await fetch(`${BASE_URL}/sessions`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
  });
  if (!res.ok) {
    throw new Error(`createSession failed: ${res.status} ${res.statusText}`);
  }
  return res.json() as Promise<SessionCreateResponse>;
}

export async function uploadH5ad(
  sessionId: string,
  file: File,
  onProgress?: (pct: number) => void
): Promise<UploadResponse> {
  const formData = new FormData();
  formData.append("file", file);

  // Signal indeterminate progress start
  onProgress?.(0);

  const res = await fetch(`${BASE_URL}/sessions/${sessionId}/upload`, {
    method: "POST",
    body: formData,
  });

  if (!res.ok) {
    throw new Error(`Upload failed: ${res.status} ${res.statusText}`);
  }

  onProgress?.(100);
  return res.json() as Promise<UploadResponse>;
}

export async function requestPlan(
  sessionId: string,
  request: PlanRequest
): Promise<Plan> {
  const res = await fetch(`${BASE_URL}/sessions/${sessionId}/plan`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(request),
  });
  if (!res.ok) {
    throw new Error(`requestPlan failed: ${res.status} ${res.statusText}`);
  }
  return res.json() as Promise<Plan>;
}

export async function executePlan(
  sessionId: string,
  plan: Plan,
  onEvent: (event: ExecuteEvent) => void
): Promise<void> {
  const res = await fetch(`${BASE_URL}/sessions/${sessionId}/execute`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ plan }),
  });

  if (!res.ok) {
    throw new Error(`executePlan failed: ${res.status} ${res.statusText}`);
  }

  if (!res.body) {
    throw new Error("No response body for SSE stream");
  }

  const reader = res.body.getReader();
  const decoder = new TextDecoder();
  let buffer = "";

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    buffer += decoder.decode(value, { stream: true });
    const lines = buffer.split("\n");
    buffer = lines.pop() ?? "";

    for (const line of lines) {
      if (line.startsWith("data: ")) {
        const raw = line.slice(6).trim();
        if (raw === "" || raw === "[DONE]") continue;
        try {
          const event = JSON.parse(raw) as ExecuteEvent;
          onEvent(event);
        } catch {
          // ignore malformed SSE lines
        }
      }
    }
  }
}

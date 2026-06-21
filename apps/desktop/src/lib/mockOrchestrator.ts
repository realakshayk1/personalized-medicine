import { http, HttpResponse } from "msw";
import type {
  SessionCreateResponse,
  UploadResponse,
  Plan,
  ExecuteEvent,
} from "@lattice/sdk/api";

const MOCK_SESSION_ID = "mock-session-001";

const MOCK_PLAN: Plan = {
  session_id: MOCK_SESSION_ID,
  workflow: "standard_scrnaseq",
  rationale:
    "Standard QC, normalization, dimensionality reduction, and clustering pipeline.",
  steps: [
    {
      primitive: "calculate_qc_metrics",
      params: { percent_top: [50, 100, 200, 500], log1p: false },
      rationale: "Compute per-cell quality metrics before filtering.",
    },
    {
      primitive: "filter_cells_basic",
      params: {
        min_genes: 200,
        max_genes: 6000,
        max_pct_mito: 15,
        min_counts: 500,
      },
      rationale: "Remove low-quality cells and likely doublets.",
    },
    {
      primitive: "normalize_total_log1p",
      params: { target_sum: null, exclude_highly_expressed: false },
      rationale: "Shifted log normalization per Theis lab best practices.",
    },
    {
      primitive: "highly_variable_genes",
      params: { n_top_genes: 2000, flavor: "seurat_v3" },
      rationale: "Select informative genes for downstream analysis.",
    },
  ],
};

const BASE = process.env["NEXT_PUBLIC_ORCHESTRATOR_URL"] ?? "http://localhost:8000";

export const handlers = [
  // POST /sessions
  http.post(`${BASE}/sessions`, () => {
    const response: SessionCreateResponse = { session_id: MOCK_SESSION_ID };
    return HttpResponse.json(response, { status: 201 });
  }),

  // POST /sessions/:id/upload
  http.post(`${BASE}/sessions/:sessionId/upload`, ({ params }) => {
    const sessionId =
      typeof params["sessionId"] === "string"
        ? params["sessionId"]
        : MOCK_SESSION_ID;
    const response: UploadResponse = {
      session_id: sessionId,
      filename: "pbmc3k.h5ad",
      sha256: "abc123def456",
      n_obs: 2700,
      n_vars: 32738,
      obs_columns: ["n_genes", "n_counts", "percent_mito", "sample", "condition"],
    };
    return HttpResponse.json(response);
  }),

  // POST /sessions/:id/plan
  http.post(`${BASE}/sessions/:sessionId/plan`, ({ params }) => {
    const sessionId =
      typeof params["sessionId"] === "string"
        ? params["sessionId"]
        : MOCK_SESSION_ID;
    const plan: Plan = { ...MOCK_PLAN, session_id: sessionId };
    return HttpResponse.json(plan);
  }),

  // POST /sessions/:id/execute  (streams SSE)
  http.post(`${BASE}/sessions/:sessionId/execute`, ({ params }) => {
    const sessionId =
      typeof params["sessionId"] === "string"
        ? params["sessionId"]
        : MOCK_SESSION_ID;

    const events: ExecuteEvent[] = [
      { type: "step_started", step_index: 0, primitive: "calculate_qc_metrics" },
      {
        type: "step_completed",
        step_index: 0,
        result: {
          primitive: "calculate_qc_metrics",
          params: { percent_top: [50, 100, 200, 500], log1p: false },
          input_hash: "abc123",
          output_hash: "def456",
          duration_sec: 1.2,
          warnings: [],
          user_explanation:
            "Computed QC metrics for 2,700 cells. Median genes per cell: 1,847.",
          methods_paragraph:
            "Quality control metrics were computed using sc.pp.calculate_qc_metrics.",
          n_obs_before: 2700,
          n_obs_after: 2700,
          n_vars_before: 32738,
          n_vars_after: 32738,
          figures: [],
        },
      },
      { type: "plan_completed", n_steps: 4 },
    ];

    const body = events
      .map((e) => `data: ${JSON.stringify({ ...e, session_id: sessionId })}\n\n`)
      .join("");

    return new HttpResponse(body, {
      headers: {
        "Content-Type": "text/event-stream",
        "Cache-Control": "no-cache",
      },
    });
  }),

  // GET /sessions/:id/provenance
  http.get(`${BASE}/sessions/:sessionId/provenance`, ({ params }) => {
    const sessionId =
      typeof params["sessionId"] === "string"
        ? params["sessionId"]
        : MOCK_SESSION_ID;
    return HttpResponse.json({
      session_id: sessionId,
      started_at: new Date().toISOString(),
      input_data: {
        filename: "pbmc3k.h5ad",
        sha256: "abc123def456",
        n_obs: 2700,
        n_vars: 32738,
      },
      steps: [],
    });
  }),
];

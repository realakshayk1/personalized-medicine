import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { FileDrop } from "../src/components/FileDrop";
import type { UploadResponse } from "@lattice/sdk/api";

describe("FileDrop", () => {
  it("renders drop zone", () => {
    render(<FileDrop sessionId="test-session" onUpload={vi.fn()} />);
    expect(screen.getByTestId("file-drop-zone")).toBeDefined();
    expect(screen.getByText(/\.h5ad/)).toBeDefined();
  });

  it("calls onUpload with response after successful drop", async () => {
    const onUpload = vi.fn();
    render(<FileDrop sessionId="mock-session-001" onUpload={onUpload} />);

    const dropZone = screen.getByTestId("file-drop-zone");
    const fakeFile = new File(["data"], "test.h5ad", {
      type: "application/octet-stream",
    });

    // Use XMLHttpRequest mock via MSW — fire a drop event
    fireEvent.drop(dropZone, {
      dataTransfer: {
        files: [fakeFile],
      },
    });

    await waitFor(
      () => {
        expect(onUpload).toHaveBeenCalledTimes(1);
        const arg = onUpload.mock.calls[0]?.[0] as UploadResponse;
        expect(arg.n_obs).toBeGreaterThan(0);
        expect(arg.n_vars).toBeGreaterThan(0);
        expect(typeof arg.session_id).toBe("string");
      },
      { timeout: 5000 }
    );
  });

  it("shows error for non-h5ad file", async () => {
    render(<FileDrop sessionId="test-session" onUpload={vi.fn()} />);
    const dropZone = screen.getByTestId("file-drop-zone");
    const badFile = new File(["data"], "test.csv", { type: "text/csv" });

    fireEvent.drop(dropZone, {
      dataTransfer: { files: [badFile] },
    });

    await waitFor(() => {
      expect(
        screen.getByText(/Only \.h5ad files are supported/)
      ).toBeDefined();
    });
  });
});

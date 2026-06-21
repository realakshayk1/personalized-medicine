import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import Home from "../src/app/page";

describe("Home page", () => {
  it("renders 4 main panels", () => {
    render(<Home />);

    expect(screen.getByTestId("panel-chat")).toBeDefined();
    expect(screen.getByTestId("figure-panel")).toBeDefined();
    expect(screen.getByTestId("panel-bottom-right")).toBeDefined();
    // chat-panel is inside panel-chat
    expect(screen.getByTestId("chat-panel")).toBeDefined();
  });

  it("renders the file drop zone", () => {
    render(<Home />);
    expect(screen.getByTestId("file-drop-zone")).toBeDefined();
  });

  it("shows Plan and Code tabs", () => {
    render(<Home />);
    expect(screen.getByText("Plan")).toBeDefined();
    expect(screen.getByText("Code")).toBeDefined();
  });
});

import { defineConfig } from "vitest/config";
import path from "path";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  test: {
    environment: "jsdom",
    setupFiles: ["./tests/setup.ts"],
    globals: true,
  },
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
      "@lattice/sdk/api": path.resolve(
        __dirname,
        "../../packages/sdk/src/api.ts"
      ),
      "@lattice/sdk": path.resolve(__dirname, "../../packages/sdk/src/index.ts"),
    },
  },
});

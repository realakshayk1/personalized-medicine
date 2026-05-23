#!/usr/bin/env bash
# Run the full Python test suite across all workspace packages.
# Invoked per-package because each package's tests/conftest.py registers under
# the same plugin name "tests.conftest", which pluggy rejects when collected
# in one session.

set -euo pipefail

echo "=== packages/primitives ==="
uv run pytest packages/primitives -q

echo "=== apps/orchestrator ==="
uv run pytest apps/orchestrator -q

echo "=== tests/ (e2e) ==="
uv run pytest tests -q

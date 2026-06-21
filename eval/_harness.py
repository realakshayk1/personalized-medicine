"""Shared eval harness utilities: result records, table printing, registry boot.

Kept separate from ``run.py`` so the suite modules (``traps.py``, ``quick.py``,
``full.py``) can import the building blocks without a circular dependency on the
argparse entry point.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

# Static run date for reproducible docs. AGENTS.md eval results must be
# reproducible, so we never call datetime.now() in a committed code path.
RUN_DATE = "2026-06-21"


@dataclass
class CaseResult:
    """Outcome of one eval case (one trap, or one quick-suite assertion)."""

    suite: str
    dataset: str
    metric: str
    value: str
    passed: bool
    detail: str = ""
    skipped: bool = False


@dataclass
class SuiteReport:
    """Aggregate outcome of a suite run."""

    suite: str
    cases: list[CaseResult] = field(default_factory=list)

    def add(self, case: CaseResult) -> None:
        self.cases.append(case)

    @property
    def n_pass(self) -> int:
        return sum(1 for c in self.cases if c.passed and not c.skipped)

    @property
    def n_fail(self) -> int:
        return sum(1 for c in self.cases if not c.passed and not c.skipped)

    @property
    def n_skip(self) -> int:
        return sum(1 for c in self.cases if c.skipped)

    @property
    def all_passed(self) -> bool:
        # Skipped cases do not fail the suite; only real failures do.
        return self.n_fail == 0


def boot_registry() -> None:
    """Import every primitive so PRIMITIVE_REGISTRY is fully populated.

    Also quiets loguru — the registry logs a DEBUG line per primitive on import
    and an INFO line per ``run_primitive`` call, which would drown the eval
    table. Eval output should be the report, not the library's chatter.
    """
    try:
        from loguru import logger

        logger.remove()
        logger.add(sys.stderr, level="WARNING")
    except Exception:  # pragma: no cover - loguru always present, defensive
        pass

    import lattice_primitives.all_primitives  # noqa: F401  (registration side effect)


def _status(case: CaseResult) -> str:
    if case.skipped:
        return "SKIP"
    return "PASS" if case.passed else "FAIL"


def print_table(report: SuiteReport) -> None:
    """Print a fixed-width results table for a suite."""
    headers = ["suite", "dataset", "metric", "value", "status"]
    rows = [[c.suite, c.dataset, c.metric, c.value, _status(c)] for c in report.cases]
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(cells: list[str]) -> str:
        return "  ".join(cell.ljust(widths[i]) for i, cell in enumerate(cells))

    sep = "  ".join("-" * w for w in widths)
    print(fmt(headers))
    print(sep)
    for row in rows:
        print(fmt(row))
    print(sep)
    print(
        f"{report.suite}: {report.n_pass} passed, {report.n_fail} failed, {report.n_skip} skipped"
    )

    # Detail lines for anything that failed or was skipped, so the reason is
    # visible without re-running with more verbosity.
    for c in report.cases:
        if (not c.passed or c.skipped) and c.detail:
            print(f"  [{_status(c)}] {c.dataset} / {c.metric}: {c.detail}")

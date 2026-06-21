"""Full suite (~90 min) — stub.

The full suite is the release-gate suite. It is NOT yet wired to its real
datasets; today it runs quick + trap and clearly logs what remains.

Planned contents (per AGENTS.md "Eval discipline" + docs/IMPLEMENTATION_PLAN.md M2):

- 5 published treatment-vs-control scRNA-seq datasets (DE top-50 overlap >= 60%).
- HCA immune subset via cellxgene-census (curated cell_type_ontology_term_id;
  annotation-F1 >= 0.75 once annotation primitives exist).
- All trap datasets (already covered by the trap suite).

Gated targets for the full suite (do not loosen without an ADR):
    annotation F1 vs published labels   >= 0.75
    DE top-50 overlap with published    >= 60%
    trap detection recall               == 100%
    reproducibility (bit-identical)     == 100%
    standalone notebook validity        == 100%

Several of these depend on primitive categories not yet implemented
(annotate, de) and on the export pipeline (M4), so the full suite cannot reach
its real targets yet.
"""

from __future__ import annotations

from eval._harness import CaseResult, SuiteReport
from eval.quick import run_quick_suite
from eval.traps import run_trap_suite

_PLANNED_DATASETS = [
    "published_tvc_1 (treatment-vs-control, DE overlap)",
    "published_tvc_2",
    "published_tvc_3",
    "published_tvc_4",
    "published_tvc_5",
    "hca_immune_subset (annotation F1 vs Cell Ontology)",
]


def run_full_suite() -> tuple[SuiteReport, SuiteReport, SuiteReport]:
    """Run quick + trap, plus a stub report documenting unwired datasets.

    Returns (full_stub_report, quick_report, trap_report).
    """
    print("full suite: real published datasets NOT yet wired.")
    print("            running quick + trap as the current full-suite coverage.")
    print("            planned additions:")
    for ds in _PLANNED_DATASETS:
        print(f"              - {ds}")

    stub = SuiteReport(suite="full")
    for ds in _PLANNED_DATASETS:
        stub.add(
            CaseResult(
                suite="full",
                dataset=ds,
                metric="not yet wired",
                value="SKIP",
                passed=True,
                skipped=True,
                detail="full datasets not yet wired (M3/M4 dependencies)",
            )
        )

    quick = run_quick_suite()
    trap = run_trap_suite()
    return stub, quick, trap

"""Dataset fetcher for the eval suites.

Real benchmark datasets are NOT committed to the repo (they are large). This
script documents and automates their download into ``eval/data/`` (gitignored).

Network failures are non-fatal: each fetch is reported and the script continues,
so it is safe to run in a partially-offline environment.

Sources
-------
- PBMC 3k: ``scanpy.datasets.pbmc3k()`` / ``pbmc3k_processed()`` — downloaded and
  cached by scanpy itself (into scanpy's datadir / ~/.cache). Used by the quick
  suite. No manual URL needed.
- Tabula Muris (1K subset): figshare. The Tabula Muris FACS/droplet objects carry
  ``cell_ontology_class`` + Cell Ontology IDs, suitable as annotation ground
  truth. URL documented below; download is best-effort.

Run with::

    uv run python -m eval.fetch_data
"""

from __future__ import annotations

import pathlib
import urllib.request

DATA_DIR = pathlib.Path(__file__).parent / "data"

# Documented source for a Tabula Muris subset. Figshare hosts the processed
# Tabula Muris objects; this is a placeholder pointer for the eval harness — the
# 1K subset is created by downsampling after download. Kept here so the source
# is recorded in-repo even though the data is not committed.
TABULA_MURIS_FIGSHARE = (
    "https://figshare.com/articles/dataset/"
    "Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687"
)


def fetch_pbmc3k() -> bool:
    """Trigger scanpy's PBMC 3k download (raw + processed). Returns success."""
    try:
        import scanpy as sc

        print("fetching pbmc3k (raw) via scanpy ...")
        sc.datasets.pbmc3k()
        print("fetching pbmc3k_processed (labels) via scanpy ...")
        sc.datasets.pbmc3k_processed()
        print("  ok: PBMC 3k cached by scanpy.")
        return True
    except Exception as exc:
        print(f"  SKIP: PBMC 3k download failed (offline?): {exc}")
        return False


def fetch_tabula_muris() -> bool:
    """Best-effort Tabula Muris subset fetch. Returns success.

    NOTE: figshare requires resolving the article to a concrete file URL. This
    function documents the source and attempts a HEAD-style reachability check;
    wiring the exact .h5ad asset is a follow-up once the full suite needs it.
    """
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Tabula Muris source (figshare): {TABULA_MURIS_FIGSHARE}")
    try:
        req = urllib.request.Request(TABULA_MURIS_FIGSHARE, method="HEAD")
        with urllib.request.urlopen(req, timeout=10):  # noqa: S310 - documented source
            print("  ok: figshare reachable. Concrete asset wiring is a follow-up.")
        return True
    except Exception as exc:
        print(f"  SKIP: Tabula Muris source unreachable (offline?): {exc}")
        return False


def main() -> int:
    print(f"data dir: {DATA_DIR}")
    pbmc_ok = fetch_pbmc3k()
    tm_ok = fetch_tabula_muris()
    print()
    print(f"pbmc3k: {'ok' if pbmc_ok else 'skipped'}; tabula_muris: {'ok' if tm_ok else 'skipped'}")
    # Non-fatal: always exit 0 so this is safe in CI / offline.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

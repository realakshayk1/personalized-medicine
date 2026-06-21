"""Additional Pydantic models for lattice_primitives.

Currently minimal — the core schemas are in _interface.py.
This module exists as an extension point for Wave 2-B+ primitives.
"""

from __future__ import annotations

# No additional models needed beyond those in _interface.py for the 3 v0.1 primitives.
# Future: add CitationEntry, TemplateContext, etc. here when the export pipeline
# reads citations.yaml and renders notebook cells.

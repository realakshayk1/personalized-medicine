"""lattice_primitives — public API.

Re-exports the public interface from _interface and registry so callers can do:

    from lattice_primitives import primitive, sanity_check, PRIMITIVE_REGISTRY, ...
"""

from lattice_primitives._interface import (
    PRIMITIVE_REGISTRY,
    AnnDataState,
    ParamSpec,
    Primitive,
    PrimitiveResult,
    Severity,
    Warning,
)
from lattice_primitives.registry import (
    PrimitiveNotFound,
    SchemaViolation,
    primitive,
    run_primitive,
    sanity_check,
)

__all__ = [
    "AnnDataState",
    "ParamSpec",
    "Primitive",
    "PrimitiveResult",
    "PrimitiveNotFound",
    "SchemaViolation",
    "Severity",
    "Warning",
    "PRIMITIVE_REGISTRY",
    "primitive",
    "run_primitive",
    "sanity_check",
]

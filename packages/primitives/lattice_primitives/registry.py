"""Primitive registry — Wave 2-A implementation.

Implements @primitive and @sanity_check decorators, SchemaViolation,
PrimitiveNotFound, and run_primitive helper used by the orchestrator.
"""

from __future__ import annotations

import functools
import hashlib
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import scipy.sparse
from anndata import AnnData
from loguru import logger

from lattice_primitives._interface import (
    PRIMITIVE_REGISTRY,
    AnnDataState,
    ParamSpec,
    Primitive,
    PrimitiveResult,
    Warning_,
)

# ---------------------------------------------------------------------------
# Exception types
# ---------------------------------------------------------------------------


class SchemaViolation(RuntimeError):
    """Raised when pre/post-condition schema validation fails."""


class PrimitiveNotFound(KeyError):
    """Raised when a primitive name is not in PRIMITIVE_REGISTRY."""


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _anndata_hash(adata: AnnData) -> str:
    """Compute a reproducible hash of AnnData.X + obs/var index."""
    h = hashlib.sha256()
    # Hash the matrix
    if scipy.sparse.issparse(adata.X):
        x_bytes = adata.X.toarray().tobytes()
    else:
        x_bytes = np.asarray(adata.X).tobytes()
    h.update(x_bytes)
    # Hash obs_names and var_names
    h.update("|".join(adata.obs_names.tolist()).encode())
    h.update("|".join(adata.var_names.tolist()).encode())
    return h.hexdigest()


def _obs_column_snapshot(adata: AnnData, columns: list[str]) -> dict[str, Any]:
    """Snapshot values of specified obs columns that exist."""
    snap: dict[str, Any] = {}
    for col in columns:
        if col in adata.obs.columns:
            snap[col] = adata.obs[col].values.copy()
    return snap


def _verify_immutable_columns(
    snapshot: dict[str, Any],
    adata_before_obs_names: list[str],
    adata_after: AnnData,
    primitive_name: str,
    row_modifying: bool,
) -> None:
    """Raise SchemaViolation if any immutable column values changed.

    For non-row-modifying primitives: requires exact shape and value match.
    For row-modifying primitives (e.g., filter_cells): checks that surviving
    cells retain their original values (subset check, not full equality).
    """
    for col, before_vals in snapshot.items():
        if col not in adata_after.obs.columns:
            raise SchemaViolation(
                f"Primitive '{primitive_name}' removed immutable obs column '{col}'"
            )
        after_vals = adata_after.obs[col].values
        before_arr = np.asarray(before_vals)
        after_arr = np.asarray(after_vals)

        if not row_modifying:
            # Shape must be identical for non-row-modifying primitives
            if before_arr.shape != after_arr.shape:
                raise SchemaViolation(
                    f"Primitive '{primitive_name}' changed shape of immutable obs column '{col}' "
                    f"({before_arr.shape} -> {after_arr.shape})"
                )
            if not np.array_equal(before_arr, after_arr):
                raise SchemaViolation(
                    f"Primitive '{primitive_name}' modified values of immutable obs column "
                    f"'{col}' (case/control swap protection triggered)"
                )
        else:
            # Row-modifying: verify surviving cells' values are unchanged
            # Build a map from obs_name -> value in before
            before_map = dict(zip(adata_before_obs_names, before_arr.tolist()))
            after_obs_names = adata_after.obs_names.tolist()
            for cell_name, after_val in zip(after_obs_names, after_arr.tolist()):
                if cell_name in before_map and before_map[cell_name] != after_val:
                    raise SchemaViolation(
                        f"Primitive '{primitive_name}' modified immutable obs column '{col}' "
                        f"for cell '{cell_name}': was '{before_map[cell_name]}', "
                        f"now '{after_val}' (case/control swap protection triggered)"
                    )


def _verify_requires(adata: AnnData, requires: AnnDataState, name: str) -> None:
    """Validate that adata satisfies the requires contract."""
    # Check obs columns
    missing_obs = [c for c in requires.obs_columns if c not in adata.obs.columns]
    if missing_obs:
        raise SchemaViolation(
            f"Primitive '{name}' requires obs columns {missing_obs} which are missing. "
            "Run the prerequisite primitive first."
        )
    # Check var columns
    missing_var = [c for c in requires.var_columns if c not in adata.var.columns]
    if missing_var:
        raise SchemaViolation(
            f"Primitive '{name}' requires var columns {missing_var} which are missing."
        )
    # Check layer_state
    if requires.layer_state is not None:
        current_state = adata.uns.get("_lattice_layer_state")
        if current_state != requires.layer_state:
            raise SchemaViolation(
                f"Primitive '{name}' requires layer_state='{requires.layer_state}' "
                f"but current state is '{current_state}'. "
                "Ensure the pipeline runs in the correct order."
            )


def _verify_produces(adata: AnnData, produces: AnnDataState, name: str) -> None:
    """Validate that adata has the columns/slots declared in produces."""
    missing_obs = [c for c in produces.obs_columns if c not in adata.obs.columns]
    if missing_obs:
        raise SchemaViolation(
            f"Primitive '{name}' declared produces.obs_columns {missing_obs} "
            "but they were not added to adata.obs."
        )
    missing_var = [c for c in produces.var_columns if c not in adata.var.columns]
    if missing_var:
        raise SchemaViolation(
            f"Primitive '{name}' declared produces.var_columns {missing_var} "
            "but they were not added to adata.var."
        )
    missing_layers = [lyr for lyr in produces.adds_layers if lyr not in adata.layers]
    if missing_layers:
        raise SchemaViolation(
            f"Primitive '{name}' declared produces.adds_layers {missing_layers} "
            "but they were not added to adata.layers."
        )
    missing_obsm = [k for k in produces.adds_obsm if k not in adata.obsm]
    if missing_obsm:
        raise SchemaViolation(
            f"Primitive '{name}' declared produces.adds_obsm {missing_obsm} "
            "but they were not added to adata.obsm."
        )
    missing_uns = [k for k in produces.adds_uns if k not in adata.uns]
    if missing_uns:
        raise SchemaViolation(
            f"Primitive '{name}' declared produces.adds_uns {missing_uns} "
            "but they were not added to adata.uns."
        )


# ---------------------------------------------------------------------------
# @primitive decorator
# ---------------------------------------------------------------------------


def primitive(
    *,
    name: str,
    category: str,
    requires: AnnDataState,
    produces: AnnDataState,
    params: dict[str, ParamSpec],
    scanpy_version: str,
    citation_key: str,
    row_modifying: bool = False,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator that registers a primitive and wraps it with schema validation.

    Parameters
    ----------
    name:
        Unique primitive name used as registry key.
    category:
        One of the allowed category literals from Primitive.
    requires:
        AnnDataState that adata must satisfy before execution.
    produces:
        AnnDataState that adata must satisfy after execution.
    params:
        Dict of parameter name -> ParamSpec.
    scanpy_version:
        Version range string for the Scanpy call inside.
    citation_key:
        Key into citations.yaml.
    row_modifying:
        If True, the primitive may change adata.n_obs (filter primitives).
        If False (default), obs_names ordering and count must be unchanged.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(func)
        def wrapper(adata: AnnData, **kwargs: Any) -> AnnData:
            logger.debug(f"Running primitive '{name}' with params {kwargs}")

            # 1. Validate requires
            _verify_requires(adata, requires, name)

            # 2. Snapshot immutable columns and shape
            imm_snapshot = _obs_column_snapshot(adata, produces.immutable_obs_columns)
            obs_names_before = adata.obs_names.tolist()
            n_obs_before = adata.n_obs

            # 3. Call the underlying function
            result_adata = func(adata, **kwargs)
            if result_adata is None:
                # in-place convention: use the original adata
                result_adata = adata

            # 4. Validate produces
            _verify_produces(result_adata, produces, name)

            # 5. Verify immutable column values are bit-identical
            if imm_snapshot:
                _verify_immutable_columns(
                    imm_snapshot,
                    obs_names_before,
                    result_adata,
                    name,
                    row_modifying,
                )

            # 6. Verify obs_names ordering if not row_modifying
            if not row_modifying:
                obs_names_after = result_adata.obs_names.tolist()
                if obs_names_after != obs_names_before:
                    raise SchemaViolation(
                        f"Primitive '{name}' is not declared row_modifying=True "
                        f"but changed obs_names (reorder or count change detected). "
                        f"Before: {n_obs_before} cells; after: {result_adata.n_obs} cells."
                    )

            # 7. Set layer_state on produces if declared
            if produces.layer_state is not None:
                result_adata.uns["_lattice_layer_state"] = produces.layer_state

            return result_adata

        # Attach primitive name for introspection
        wrapper.__primitive_name__ = name  # type: ignore[attr-defined]

        # Resolve template path
        import pathlib

        pkg_dir = pathlib.Path(__file__).parent
        template_candidate = pkg_dir / "templates" / f"{name}.j2"
        template_path = str(template_candidate) if template_candidate.exists() else None

        # Build the Primitive registry entry
        entry = Primitive(
            name=name,
            category=category,  # type: ignore[arg-type]
            requires=requires,
            produces=produces,
            params=params,
            scanpy_version=scanpy_version,
            citation_key=citation_key,
            func=wrapper,
            sanity_check_func=None,
            template_path=template_path,
        )
        PRIMITIVE_REGISTRY[name] = entry
        logger.debug(f"Registered primitive '{name}'")

        return wrapper

    return decorator


# ---------------------------------------------------------------------------
# @sanity_check decorator
# ---------------------------------------------------------------------------


def sanity_check(
    primitive_name: str,
) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator that attaches a sanity-check function to a registered primitive.

    The decorated function must have the signature:
        (adata_before: AnnData, adata_after: AnnData, params: dict) -> list[Warning_]

    Registration is deferred — the primitive does not need to be registered
    before the sanity_check decorator runs, as long as both are imported before
    run_primitive is called.
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        if primitive_name in PRIMITIVE_REGISTRY:
            # Rebuild the frozen Primitive with the sanity_check_func attached
            entry = PRIMITIVE_REGISTRY[primitive_name]
            PRIMITIVE_REGISTRY[primitive_name] = entry.model_copy(
                update={"sanity_check_func": func}
            )
        else:
            # Store for deferred attachment — will be attached on first run_primitive
            _PENDING_SANITY_CHECKS[primitive_name] = func
        return func

    return decorator


# Deferred sanity-check registry for primitives not yet imported
_PENDING_SANITY_CHECKS: dict[str, Callable[..., Any]] = {}


def _flush_pending_sanity_checks() -> None:
    """Attach any pending sanity checks to now-registered primitives."""
    for pname, func in list(_PENDING_SANITY_CHECKS.items()):
        if pname in PRIMITIVE_REGISTRY:
            entry = PRIMITIVE_REGISTRY[pname]
            PRIMITIVE_REGISTRY[pname] = entry.model_copy(
                update={"sanity_check_func": func}
            )
            del _PENDING_SANITY_CHECKS[pname]


# ---------------------------------------------------------------------------
# run_primitive helper
# ---------------------------------------------------------------------------


def run_primitive(
    name: str,
    adata: AnnData,
    **params: Any,
) -> tuple[AnnData, list[Warning_], PrimitiveResult]:
    """Execute a named primitive and return (adata, warnings, result).

    This is the entry point used by the orchestrator. It:
    1. Looks up the primitive in PRIMITIVE_REGISTRY.
    2. Flushes any pending sanity-check registrations.
    3. Computes input hash.
    4. Calls the registered (wrapped) function.
    5. Calls the sanity_check_func if present.
    6. Computes output hash.
    7. Renders explanation templates if available.
    8. Returns a structured PrimitiveResult.
    """
    _flush_pending_sanity_checks()

    if name not in PRIMITIVE_REGISTRY:
        raise PrimitiveNotFound(
            f"Primitive '{name}' not found in PRIMITIVE_REGISTRY. "
            f"Available: {list(PRIMITIVE_REGISTRY.keys())}"
        )

    entry = PRIMITIVE_REGISTRY[name]

    n_obs_before = adata.n_obs
    n_vars_before = adata.n_vars
    input_hash = _anndata_hash(adata)

    start = time.perf_counter()
    adata_out = entry.func(adata, **params)
    duration_sec = time.perf_counter() - start

    output_hash = _anndata_hash(adata_out)
    n_obs_after = adata_out.n_obs
    n_vars_after = adata_out.n_vars

    # Run sanity checks
    # Pass n_obs_before/n_vars_before in params so checks can compare against
    # pre-execution counts even when the primitive mutates adata in-place
    # (in which case adata_before IS adata_after after the call).
    sanity_params = {
        **params,
        "_n_obs_before": n_obs_before,
        "_n_vars_before": n_vars_before,
    }
    warnings: list[Warning_] = []
    if entry.sanity_check_func is not None:
        try:
            warnings = entry.sanity_check_func(adata, adata_out, sanity_params)
        except Exception as exc:
            logger.warning(f"Sanity check for '{name}' raised an exception: {exc}")
            warnings = [
                Warning_(
                    severity="warn",
                    message=f"Sanity check errored: {exc}",
                    primitive=name,
                )
            ]

    # Render Jinja templates if available
    user_explanation = ""
    methods_paragraph = ""
    if entry.template_path is not None:
        try:
            user_explanation, methods_paragraph = _render_templates(
                entry.template_path, name, params, adata_out
            )
        except Exception as exc:
            logger.warning(f"Template rendering for '{name}' failed: {exc}")

    result = PrimitiveResult(
        primitive=name,
        params=params,
        input_hash=input_hash,
        output_hash=output_hash,
        duration_sec=duration_sec,
        warnings=warnings,
        user_explanation=user_explanation,
        methods_paragraph=methods_paragraph,
        n_obs_before=n_obs_before,
        n_obs_after=n_obs_after,
        n_vars_before=n_vars_before,
        n_vars_after=n_vars_after,
    )

    logger.info(
        f"Primitive '{name}' completed in {duration_sec:.3f}s "
        f"({n_obs_before}->{n_obs_after} cells, {n_vars_before}->{n_vars_after} genes). "
        f"{len(warnings)} warning(s)."
    )

    return adata_out, warnings, result


def _render_templates(
    template_path: str,
    name: str,
    params: dict[str, Any],
    adata: AnnData,
) -> tuple[str, str]:
    """Render Jinja template and return (user_block, methods_block)."""
    import pathlib

    from jinja2 import Environment, FileSystemLoader

    tpath = pathlib.Path(template_path)
    env = Environment(loader=FileSystemLoader(str(tpath.parent)), autoescape=False)
    tmpl = env.get_template(tpath.name)

    context = {
        "name": name,
        "params": params,
        "adata": adata,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
    }

    # Extract each block separately for reliability
    user_block = _extract_block(tmpl, "user", context)
    methods_block = _extract_block(tmpl, "methods", context)

    return user_block.strip(), methods_block.strip()


def _extract_block(tmpl: Any, block_name: str, context: dict[str, Any]) -> str:
    """Render only one named block from a Jinja2 template."""
    try:
        rendered_blocks = tmpl.blocks.get(block_name)
        if rendered_blocks:
            ctx = tmpl.new_context(context)
            return "".join(rendered_blocks[0](ctx))
    except Exception:
        pass
    return ""

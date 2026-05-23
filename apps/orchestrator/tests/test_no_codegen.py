"""Test that the orchestrator package contains no free-form code generation paths.

Per ADR-001, the orchestrator MUST NOT contain:
- exec(...) calls on user-derived strings
- eval(...) calls on user-derived strings
- subprocess.run / subprocess.call / subprocess.Popen with non-literal first args
  that could execute user-derived code

This test uses the AST to parse every .py file under lattice/ and inspects
Call nodes — it only rejects actual Call nodes in code, not comments or strings.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

LATTICE_DIR = pathlib.Path(__file__).parent.parent / "lattice"

# Names that are forbidden when called
FORBIDDEN_CALLABLES = {"exec", "eval"}

# subprocess functions that could be exploited
SUBPROCESS_FUNCS = {"run", "call", "check_call", "check_output", "Popen"}


def _is_string_literal(node: ast.expr) -> bool:
    """Return True if node is a string or bytes literal."""
    return isinstance(node, ast.Constant) and isinstance(node.value, (str, bytes))


def _is_list_of_literals(node: ast.expr) -> bool:
    """Return True if node is a list/tuple of string literals (safe subprocess args)."""
    if isinstance(node, (ast.List, ast.Tuple)):
        return all(_is_string_literal(elt) for elt in node.elts)
    return False


class ForbiddenCallVisitor(ast.NodeVisitor):
    """Finds forbidden call patterns in an AST."""

    def __init__(self) -> None:
        self.violations: list[str] = []

    def visit_Call(self, node: ast.Call) -> None:  # noqa: N802
        # Check for exec(...) or eval(...) as direct names
        if isinstance(node.func, ast.Name) and node.func.id in FORBIDDEN_CALLABLES:
            self.violations.append(
                f"Line {node.lineno}: direct call to {node.func.id}()"
            )

        # Check for builtins.exec / builtins.eval
        if (
            isinstance(node.func, ast.Attribute)
            and node.func.attr in FORBIDDEN_CALLABLES
        ):
            self.violations.append(
                f"Line {node.lineno}: attribute call {node.func.attr}()"
            )

        # Check for subprocess.run/call/etc. with non-literal first arg
        if isinstance(node.func, ast.Attribute) and node.func.attr in SUBPROCESS_FUNCS:
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "subprocess":
                # Check the first positional argument
                if node.args:
                    first_arg = node.args[0]
                    if not (_is_string_literal(first_arg) or _is_list_of_literals(first_arg)):
                        self.violations.append(
                            f"Line {node.lineno}: subprocess.{node.func.attr}() "
                            f"with non-literal first arg (potential code injection)"
                        )
                elif node.keywords:
                    # args= keyword
                    for kw in node.keywords:
                        if kw.arg == "args" and not (
                            _is_string_literal(kw.value) or _is_list_of_literals(kw.value)
                        ):
                            self.violations.append(
                                f"Line {node.lineno}: subprocess.{node.func.attr}(args=...) "
                                f"with non-literal args (potential code injection)"
                            )

        self.generic_visit(node)


def _check_file(path: pathlib.Path) -> list[str]:
    """Parse a Python file and return a list of violations."""
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as exc:
        return [f"SyntaxError: {exc}"]

    visitor = ForbiddenCallVisitor()
    visitor.visit(tree)
    return [f"{path}: {v}" for v in visitor.violations]


def test_no_codegen_patterns() -> None:
    """Grep the lattice/ package AST for forbidden code-gen patterns."""
    all_violations: list[str] = []

    py_files = list(LATTICE_DIR.rglob("*.py"))
    assert len(py_files) > 0, f"No Python files found under {LATTICE_DIR}"

    for py_file in py_files:
        violations = _check_file(py_file)
        all_violations.extend(violations)

    if all_violations:
        violation_str = "\n".join(all_violations)
        pytest.fail(
            f"Found {len(all_violations)} forbidden code-gen pattern(s) in lattice/ "
            f"(ADR-001 violation):\n{violation_str}"
        )

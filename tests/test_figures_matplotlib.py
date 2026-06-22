from __future__ import annotations

import ast
from pathlib import Path


FIGURE_CORE = Path(__file__).resolve().parents[1] / "astromol" / "figures" / "_core.py"


def _keyword_value(call: ast.Call, name: str) -> ast.expr | None:
    for keyword in call.keywords:
        if keyword.arg == name:
            return keyword.value
    return None


def test_named_pyplot_figures_clear_existing_figures() -> None:
    tree = ast.parse(FIGURE_CORE.read_text())
    missing_clear: list[int] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        function = node.func
        if not (
            isinstance(function, ast.Attribute)
            and isinstance(function.value, ast.Name)
            and function.value.id == "plt"
            and function.attr in {"figure", "subplots"}
        ):
            continue
        if _keyword_value(node, "num") is None:
            continue
        clear_value = _keyword_value(node, "clear")
        if not (
            isinstance(clear_value, ast.Constant)
            and clear_value.value is True
        ):
            missing_clear.append(node.lineno)

    assert missing_clear == []

# -*- coding: utf-8 -*-
"""Import helpers for bundled packages without mutating sys.path."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def ensure_bundled_mhm_tools():
    """Register the bundled ``mhm_tools`` package as a top-level package."""
    package_root = Path(__file__).resolve().parents[1] / "mhm_tools"
    init_file = package_root / "__init__.py"
    if not init_file.exists():
        return None

    existing = sys.modules.get("mhm_tools")
    if existing is not None:
        existing_file = getattr(existing, "__file__", "")
        if existing_file and Path(existing_file).resolve() == init_file.resolve():
            return existing
        if getattr(existing, "__path__", None):
            return existing

    spec = importlib.util.spec_from_file_location(
        "mhm_tools",
        str(init_file),
        submodule_search_locations=[str(package_root)],
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load bundled mhm_tools from {init_file}")

    module = importlib.util.module_from_spec(spec)
    sys.modules["mhm_tools"] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        if sys.modules.get("mhm_tools") is module:
            sys.modules.pop("mhm_tools", None)
        raise
    return module


__all__ = ["ensure_bundled_mhm_tools"]

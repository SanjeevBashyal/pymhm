# -*- coding: utf-8 -*-
"""Running QGIS Processing algorithms.

The one place that imports `processing`. That import only resolves inside a
real QGIS session -- the standalone shim does not provide a Processing backend
-- so keeping it here means a module that merely *sequences* algorithms stays
importable outside QGIS.

`run` never raises: a failed algorithm returns None, having reported through
the callbacks, because every caller already branches on a falsy result.
"""
from __future__ import annotations

from ...core.utils.report import CRITICAL


def run(name: str, params: dict, *, log=None, report=None):
    """Run one Processing algorithm, returning None when it fails."""
    if log is not None:
        log(f"Running algorithm: {name}...")
    try:
        import processing
    except ImportError as error:
        message = f"QGIS Processing is unavailable: {error}"
        if log is not None:
            log(f"ERROR: {message}")
        if report is not None:
            report(CRITICAL, "Processing Error", message)
        return None

    try:
        result = processing.run(name, params)
    except Exception as error:
        if log is not None:
            log(f"ERROR: Algorithm '{name}' failed. Details: {error}")
        if report is not None:
            report(CRITICAL, "Processing Error",
                   f"Algorithm '{name}' failed.\nCheck the log for details.")
        return None

    if log is not None:
        log(f"Algorithm '{name}' finished successfully.")
    return result


__all__ = ["run"]

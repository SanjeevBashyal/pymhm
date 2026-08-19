# -*- coding: utf-8 -*-
"""Reuse checks for meteorology forcing already recorded in the project state.

A prepared file may only be reused when the processing state records it, the
file and its header are still on disk, and that header still matches the L2
grid the current run targets. A changed L2 grid always forces a rebuild.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Iterable, Mapping

from .l2_grid import assert_header_file_matches
from .paths import expected_meteo_outputs
from ..project_layout import workspace_folder


BASE_VARIABLES = ("pre", "tavg", "tmin", "tmax")


def required_meteo_variables(include_pet: bool) -> tuple[str, ...]:
    """Return the forcing variables one run has to produce."""
    if include_pet:
        return BASE_VARIABLES + ("pet",)
    return BASE_VARIABLES


def output_state_key(project_folder, path) -> str:
    """Return the processing-state key recorded for a meteorology output."""
    try:
        return os.path.relpath(
            str(path), workspace_folder(project_folder)
        ).replace("\\", "/")
    except ValueError:
        return str(Path(path).resolve()).replace("\\", "/")


def stale_meteo_variables(
        project_folder,
        recorded_outputs: Mapping[str, Any],
        header: Mapping[str, Any],
        variables: Iterable[str]) -> dict[str, str]:
    """Return the variables that still need preparing, mapped to the reason."""
    expected = expected_meteo_outputs(project_folder)
    stale: dict[str, str] = {}
    for variable in variables:
        paths = expected.get(variable)
        if paths is None:
            stale[variable] = "unknown meteorology variable"
            continue
        for path in (paths["netcdf"], paths["header"]):
            key = output_state_key(project_folder, path)
            entry = recorded_outputs.get(key)
            if not isinstance(entry, Mapping) or not entry.get("exists"):
                stale[variable] = f"{path.name} is not recorded as prepared"
                break
            if not Path(path).is_file():
                stale[variable] = f"{path.name} is missing on disk"
                break
        else:
            try:
                assert_header_file_matches(paths["header"], header)
            except (OSError, ValueError) as error:
                stale[variable] = f"{error}"
    return stale


__all__ = [
    "output_state_key",
    "required_meteo_variables",
    "stale_meteo_variables",
]

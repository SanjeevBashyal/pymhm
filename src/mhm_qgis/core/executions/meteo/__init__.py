# -*- coding: utf-8 -*-
"""Preparing mHM meteorology forcing.

Qt-free: the caller passes `log=` for the running commentary and `report=` for
anything a user has to be shown, so the same functions serve the dialog, the
worker thread and the tests.
"""
from __future__ import annotations

import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from ...meteorology.forcing import (
    MeteoFolderSpec,
    SpatialMetadata,
    process_meteo_inputs,
)
from ...grid import TargetGrid
from ...meteorology.mask import write_meteo_mask
from ...handlers.store.layout import (
    expected_meteo_outputs,
    is_v6,
    meteo_mask_path,
    meteo_output_root,
)
from ...meteorology.reuse import required_meteo_variables, stale_meteo_variables
from ...handlers.state import processing
from ...utils.report import CRITICAL, INFORMATION


@dataclass(frozen=True)
class MeteorologyRun:
    """Immutable inputs prepared on the UI thread for one meteo run."""

    project_folder: Path
    precipitation: MeteoFolderSpec
    temperature: MeteoFolderSpec
    pet: MeteoFolderSpec | None
    target_grid: TargetGrid
    grid_metadata: dict
    inspections: Mapping[str, SpatialMetadata] | None = None


def _say(log, message):
    if log is not None:
        log(message)


def existing_outputs(project_folder, state, *, log=None) -> list[str]:
    """Record the forcing files already on disk and return their variable names."""
    if not project_folder:
        return []

    found = []
    for variable, paths in expected_meteo_outputs(project_folder).items():
        header = paths["header"]
        if not paths["netcdf"].exists():
            continue
        if header is not None and not header.exists():
            continue
        state.mark_prepared(paths["netcdf"], paths["netcdf"].name)
        if header is not None:
            state.mark_prepared(header, header.name)
        found.append(variable)

    if found:
        _say(log, "Existing meteorology forcing files found: " + ", ".join(found))
    return found


def reuse_prepared(run: MeteorologyRun, state, *, log=None) -> bool:
    """Return True when recorded forcing already matches the target L2 grid."""
    variables = required_meteo_variables(run.pet is not None)
    stale = stale_meteo_variables(
        run.project_folder,
        state.load().get("outputs", {}),
        run.target_grid.header,
        variables,
    )
    if stale:
        for variable, reason in sorted(stale.items()):
            _say(log, f"{variable}: needs preparing ({reason}).")
        return False

    _say(log,
         "Meteorology forcing is already prepared on this L2 grid for "
         f"{', '.join(variables)}. Skipping meteorology preparation.")
    _say(log,
         "Change the L2 grid or delete the files under "
         f"{meteo_output_root(run.project_folder)} to force a rebuild.")
    processing.save_grid_metadata(run.project_folder, run.grid_metadata)
    return True


def prepare_forcing(run: MeteorologyRun, state, *, log=None, report=None) -> bool:
    """Prepare every selected meteorology forcing file. Returns success."""
    _say(log, "\n--- Preparing meteorology forcing files ---")
    _say(log, f"Meteorology output folder: {meteo_output_root(run.project_folder)}")
    if reuse_prepared(run, state, log=log):
        return True

    try:
        flat = is_v6(run.project_folder)
        result = process_meteo_inputs(
            precipitation=run.precipitation,
            temperature=run.temperature,
            pet=run.pet,
            output_root=meteo_output_root(run.project_folder),
            target_grid=run.target_grid,
            flat_layout=flat,
            log=log,
            inspections=run.inspections,
        )
        if flat:
            # v6 names a meteorology mask in config_input.
            mask = write_meteo_mask(meteo_mask_path(run.project_folder), run.target_grid)
            state.mark_prepared(mask, mask.name)
            _say(log, f"Meteorology mask written: {mask}")
    except Exception as error:
        _say(log, f"ERROR: Meteorology processing failed: {error}\n{traceback.format_exc()}")
        if report is not None:
            report(CRITICAL, "Meteorology Error",
                   f"Meteorology processing failed.\n{error}")
        return False

    for variable, output_path in result.outputs.items():
        state.mark_prepared(output_path, Path(output_path).name)
        header_path = result.headers.get(variable)
        if header_path:
            state.mark_prepared(header_path, Path(header_path).name)
    _say(log, "Meteorology forcing preparation completed.")
    _say(log, f"Prepared variables: {', '.join(sorted(result.outputs))}")
    processing.save_grid_metadata(run.project_folder, run.grid_metadata)
    if report is not None:
        report(INFORMATION, "Success",
               "Meteorology forcing files prepared successfully.")
    return True


__all__ = [
    "MeteorologyRun",
    "existing_outputs",
    "prepare_forcing",
    "reuse_prepared",
]

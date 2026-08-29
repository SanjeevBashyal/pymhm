# -*- coding: utf-8 -*-
"""Meteorology run input collection and validation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .domain import bounds_from_dialog
from .exceptions import MeteorologyInputError
from .paths import meteo_output_root


@dataclass(frozen=True)
class MeteorologyRunInputs:
    """Validated inputs needed to prepare mHM meteorology forcing files."""

    project_folder: Path
    nc_folder: Path
    output_root: Path
    crop_bounds: tuple[float, float, float, float] | None


def collect_run_inputs(dialog: object) -> MeteorologyRunInputs:
    """Read the dialog and return validated meteorology processing inputs."""
    if not dialog.project_folder:
        raise MeteorologyInputError(
            "Missing Input",
            "Please select a project folder before preparing meteorology data.",
            severity="critical",
        )

    meteo_folder_text = dialog.lineEdit_meteo_folder.text().strip()
    if not meteo_folder_text:
        raise MeteorologyInputError(
            "Missing Input",
            "Please select the folder containing ERA5-Land NetCDF files.",
        )

    nc_folder = Path(meteo_folder_text)
    if not nc_folder.exists() or not nc_folder.is_dir():
        raise MeteorologyInputError(
            "Invalid Input",
            f"Meteorology data folder does not exist:\n{nc_folder}",
            severity="critical",
        )

    bounds = bounds_from_dialog(dialog)
    crop_bounds = None if bounds is None else bounds.as_tuple()
    project_folder = Path(dialog.project_folder)
    return MeteorologyRunInputs(
        project_folder=project_folder,
        nc_folder=nc_folder,
        output_root=meteo_output_root(project_folder),
        crop_bounds=crop_bounds,
    )

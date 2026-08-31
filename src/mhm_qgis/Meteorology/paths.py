# -*- coding: utf-8 -*-
"""Meteorology output path helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Union
from ..core.handlers.store.paths import meteo_folder
from ..core.handlers.store.layout import is_v6


METEO_VARIABLES = ("pre", "tavg", "tmin", "tmax", "pet")
PathInput = Union[str, Path]


def meteo_output_root(project_folder: PathInput) -> Path:
    """Return the project-local mHM meteorology forcing root."""
    return Path(meteo_folder(project_folder))


def meteo_mask_path(project_folder: PathInput) -> Path:
    """Return the v6 meteorology mask path."""
    return meteo_output_root(project_folder) / "mask.nc"


def expected_meteo_outputs(project_folder: PathInput) -> dict[str, dict[str, Path]]:
    """Return expected final meteo files keyed by variable.

    v6 keeps one flat NetCDF per variable in `data/master/meteo` and has no
    accompanying header; v5.13 keeps a per-variable folder with `header.txt`.
    """
    root = meteo_output_root(project_folder)
    if is_v6(project_folder):
        return {
            variable: {"netcdf": root / f"{variable}.nc", "header": None}
            for variable in METEO_VARIABLES
        }
    return {
        variable: {
            "netcdf": root / variable / f"{variable}.nc",
            "header": root / variable / "header.txt",
        }
        for variable in METEO_VARIABLES
    }

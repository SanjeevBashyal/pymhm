# -*- coding: utf-8 -*-
"""Meteorology output path helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Union

from ..project_layout import meteo_folder


METEO_VARIABLES = ("pre", "tavg", "tmin", "tmax")
PathInput = Union[str, Path]


def meteo_output_root(project_folder: PathInput) -> Path:
    """Return the project-local mHM meteorology forcing root."""
    return Path(meteo_folder(project_folder))


def expected_meteo_outputs(project_folder: PathInput) -> dict[str, dict[str, Path]]:
    """Return expected final meteo NetCDF and header files keyed by variable."""
    root = meteo_output_root(project_folder)
    return {
        variable: {
            "netcdf": root / variable / f"{variable}.nc",
            "header": root / variable / "header.txt",
        }
        for variable in METEO_VARIABLES
    }

# -*- coding: utf-8 -*-
"""Resolve prepared meteorology forcing for the display selector.

Kept free of QGIS so the band arithmetic stays testable; the canvas half lives
in :mod:`mhm_qgis.qgis_bridge.display.meteo`.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta
from pathlib import Path

from ..handlers.store.layout import METEO_VARIABLES, expected_meteo_outputs

#: The forcing time axis is stored as days since this epoch.
TIME_EPOCH = date(1900, 1, 1)

#: `pet` is optional, so the display order lists it last.
DISPLAY_VARIABLES = ("pre", "tavg", "tmin", "tmax", "pet")

VARIABLE_LABELS = {
    "pre": "Precipitation",
    "tavg": "Mean temperature",
    "tmin": "Minimum temperature",
    "tmax": "Maximum temperature",
    "pet": "Potential evapotranspiration",
}


@dataclass(frozen=True)
class MeteoOutput:
    """One timestep of one prepared forcing variable."""

    path: Path
    variable: str
    band: int
    when: date
    name: str


def meteo_output_path(project_folder, variable: str) -> Path | None:
    """Return the prepared NetCDF for one variable, when it exists on disk."""
    paths = expected_meteo_outputs(project_folder).get(variable)
    if paths is None:
        return None
    netcdf = Path(paths["netcdf"])
    return netcdf if netcdf.is_file() else None


def available_meteo_variables(project_folder) -> list[str]:
    """Return the forcing variables actually prepared for this project.

    `expected_meteo_outputs` names every variable including `pet`, which is only
    written when the run asked for it, so each candidate is checked on disk.
    """
    return [
        variable
        for variable in DISPLAY_VARIABLES
        if variable in METEO_VARIABLES
        and meteo_output_path(project_folder, variable) is not None
    ]


def _time_axis_bounds(path) -> tuple[int, int] | None:
    """Return the first and last raw day offsets, without reading the axis."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        return None
    try:
        with Dataset(str(path)) as dataset:
            axis = dataset.variables.get("time")
            if axis is None or axis.size == 0:
                return None
            return int(axis[0]), int(axis[-1])
    except OSError:
        return None


def meteo_time_range(path) -> tuple[date, date] | None:
    """Return the first and last date covered by a prepared forcing file."""
    bounds = _time_axis_bounds(path)
    if bounds is None:
        return None
    first, last = bounds
    return TIME_EPOCH + timedelta(days=first), TIME_EPOCH + timedelta(days=last)


def project_time_range(project_folder) -> tuple[date, date] | None:
    """Return the forcing period, measured on the first prepared variable."""
    for variable in available_meteo_variables(project_folder):
        span = meteo_time_range(meteo_output_path(project_folder, variable))
        if span is not None:
            return span
    return None


def step_count(span: tuple[date, date]) -> int:
    """Return the number of daily steps a forcing period covers."""
    first, last = span
    return (last - first).days + 1


def date_for_step(span: tuple[date, date], step: int) -> date:
    """Return the date at one zero-based step, clamped to the period."""
    first, last = span
    return min(max(first + timedelta(days=int(step)), first), last)


def step_for_date(span: tuple[date, date], when) -> int:
    """Return the zero-based step for one date, clamped to the period."""
    first, last = span
    when = when.date() if isinstance(when, datetime) else when
    return max(0, min((when - first).days, (last - first).days))


def resolve_meteo_output(project_folder, variable: str, when) -> MeteoOutput | None:
    """Return the raster band for one variable at one date, if it is covered.

    The axis is daily and sorted, so the band follows from arithmetic rather
    than a decode of the whole series.
    """
    path = meteo_output_path(project_folder, variable)
    if path is None:
        return None
    bounds = _time_axis_bounds(path)
    if bounds is None:
        return None
    first, last = bounds
    when = when.date() if isinstance(when, datetime) else when
    offset = (when - TIME_EPOCH).days
    if not first <= offset <= last:
        return None
    label = VARIABLE_LABELS.get(variable, variable)
    return MeteoOutput(
        path=path,
        variable=variable,
        band=offset - first + 1,
        when=when,
        name=f"{label} ({variable})",
    )


__all__ = [
    "DISPLAY_VARIABLES",
    "MeteoOutput",
    "VARIABLE_LABELS",
    "available_meteo_variables",
    "date_for_step",
    "meteo_output_path",
    "meteo_time_range",
    "project_time_range",
    "resolve_meteo_output",
    "step_count",
    "step_for_date",
]

# -*- coding: utf-8 -*-
"""ERA5-Land to mHM forcing API independent of the QGIS UI."""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from ..logging import capture_messages
from .era5_mhm import (
    MissingDependencyError,
    inspect_era5_folder,
    process_era5_to_mhm as _process_era5_to_mhm,
)
from .era5_mhm.types import MeteoForcingResult


def process_era5_to_mhm(
        nc_folder: str | Path,
        output_root: str | Path,
        bounds: tuple[float, float, float, float] | None = None,
        target_lat=None,
        target_lon=None,
        target_header: dict | None = None,
        skip_existing: bool = True,
        log: Callable[[str], None] | None = None) -> MeteoForcingResult:
    """Prepare mHM meteorology forcing files from ERA5-Land NetCDF files."""
    with capture_messages(log):
        return _process_era5_to_mhm(
            nc_folder=nc_folder,
            output_root=output_root,
            bounds=bounds,
            target_lat=target_lat,
            target_lon=target_lon,
            target_header=target_header,
            skip_existing=skip_existing,
            log=log,
        )


__all__ = [
    "MeteoForcingResult",
    "MissingDependencyError",
    "inspect_era5_folder",
    "process_era5_to_mhm",
]

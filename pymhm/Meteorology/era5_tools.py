# -*- coding: utf-8 -*-
"""Lazy import adapter for ERA5-Land processing tools."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Any

from .exceptions import MeteorologyToolImportError


@dataclass(frozen=True)
class ERA5MhmTools:
    """Callable ERA5-Land tools used by the QGIS meteorology workflow."""

    inspect_folder: Callable[[Path], dict[str, Any]]
    process_to_mhm: Callable[..., Any]
    missing_dependency_error: type[Exception]


def load_era5_mhm_tools() -> ERA5MhmTools:
    """Load ERA5-Land processing functions without importing them at plugin load."""
    try:
        from .ERA5Land import inspect_era5_folder, process_era5_to_mhm
        from .ERA5Land.mhm_forcing import MissingDependencyError
    except Exception as e:
        raise MeteorologyToolImportError(
            f"Could not import ERA5-Land processing tools:\n{e}") from e

    return ERA5MhmTools(
        inspect_folder=inspect_era5_folder,
        process_to_mhm=process_era5_to_mhm,
        missing_dependency_error=MissingDependencyError,
    )

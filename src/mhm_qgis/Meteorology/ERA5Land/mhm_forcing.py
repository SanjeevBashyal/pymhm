"""Backward-compatible facade for ERA5-Land to mHM forcing conversion."""
from __future__ import annotations

from .mhm import (
    FORCING_SPECS,
    ForcingSpec,
    MeteoForcingResult,
    MissingDependencyError,
    inspect_era5_folder,
    process_era5_to_mhm,
)

__all__ = [
    "FORCING_SPECS",
    "ForcingSpec",
    "MeteoForcingResult",
    "MissingDependencyError",
    "inspect_era5_folder",
    "process_era5_to_mhm",
]

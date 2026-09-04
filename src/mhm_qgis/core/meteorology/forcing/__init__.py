# -*- coding: utf-8 -*-
"""Inspect inputs and prepare model-ready meteorology forcing."""

from .api import (
    MeteoForcingResult,
    inspect_meteo_folder,
    inspect_meteo_inputs,
    process_meteo_inputs,
    resolution_in_crs,
)
from .cache import (
    inspect_meteo_folder_cached,
    inspect_meteo_inputs_cached,
    inspection_fingerprint,
)
from .types import (
    ERA5LAND,
    MHM_READY,
    MeteoFolderSpec,
    SpatialMetadata,
    SpatialResolution,
    TargetGrid,
    normalize_kind,
    normalize_source,
)

__all__ = [
    "ERA5LAND",
    "MHM_READY",
    "MeteoFolderSpec",
    "MeteoForcingResult",
    "SpatialMetadata",
    "SpatialResolution",
    "TargetGrid",
    "inspect_meteo_folder",
    "inspect_meteo_folder_cached",
    "inspect_meteo_inputs",
    "inspect_meteo_inputs_cached",
    "inspection_fingerprint",
    "normalize_kind",
    "normalize_source",
    "process_meteo_inputs",
    "resolution_in_crs",
]

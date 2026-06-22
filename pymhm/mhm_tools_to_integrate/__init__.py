# -*- coding: utf-8 -*-
"""UI-free core integration package for pymhm workflows."""
from __future__ import annotations

from ._bundled import ensure_bundled_mhm_tools

from .data_processing import (
    MissingDependencyError,
    MorphologyAsciiLayer,
    MorphologyAsciiResult,
    RegionDataHandler,
    prepare_morphology_ascii_files,
    inspect_era5_folder,
    process_era5_to_mhm,
    validate_grid_headers,
)
from .logging import capture_messages
from .setup_creation import (
    Domain,
    Outlet,
    Region,
    TerrainProducts,
    Watershed,
    create_terrain_products,
)
from .setup_creation.latlon import create_latlon_file

__all__ = [
    "Domain",
    "MissingDependencyError",
    "MorphologyAsciiLayer",
    "MorphologyAsciiResult",
    "Outlet",
    "Region",
    "RegionDataHandler",
    "TerrainProducts",
    "Watershed",
    "capture_messages",
    "create_latlon_file",
    "create_terrain_products",
    "ensure_bundled_mhm_tools",
    "inspect_era5_folder",
    "prepare_morphology_ascii_files",
    "process_era5_to_mhm",
    "validate_grid_headers",
]

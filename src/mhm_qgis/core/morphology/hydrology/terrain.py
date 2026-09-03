# -*- coding: utf-8 -*-
"""Terrain derivatives: aspect and slope.

`process_aspect` and `process_slope` were two near-identical methods differing
only in output name, algorithm and slope's scale factor. They are one function
over a product description here.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable

from ...handlers.store.paths import geometry_folder
from ...utils.report import WARNING

#: gdaldem wants a degrees-to-metres factor for a geographic DEM; a projected
#: one is already in metres.
GEOGRAPHIC_SLOPE_SCALE = 111200.0


def aspect_params(dem_path, output_path):
    """Return GDAL aspect parameters."""
    return {
        "INPUT": str(dem_path),
        "BAND": 1,
        "TRIG_ANGLE": False,
        "ZERO_FLAT": False,
        "COMPUTE_EDGES": False,
        "ZEVENBERGEN": False,
        "OPTIONS": None,
        "EXTRA": "",
        "OUTPUT": str(output_path),
    }


def slope_params(dem_path, output_path, scale=1.0):
    """Return GDAL slope parameters."""
    return {
        "INPUT": str(dem_path),
        "BAND": 1,
        "SCALE": scale,
        "AS_PERCENT": True,
        "COMPUTE_EDGES": False,
        "ZEVENBERGEN": False,
        "OPTIONS": None,
        "EXTRA": "",
        "OUTPUT": str(output_path),
    }


@dataclass(frozen=True)
class TerrainProduct:
    """One gdaldem derivative: where it lands and how it is produced."""

    field: str          # the session attribute holding its path
    filename: str
    layer_name: str
    algorithm: str
    label: str


ASPECT = TerrainProduct("aspect_path", "1_dem_aspect.tif", "1_DEM_Aspect",
                        "gdal:aspect", "Aspect")
SLOPE = TerrainProduct("slope_path", "1_dem_slope.tif", "1_DEM_Slope",
                       "gdal:slope", "Slope")


def slope_scale(crs) -> float:
    """Match gdaldem's scale choice: geographic DEMs need degrees to metres."""
    if crs is not None and crs.isValid() and crs.isGeographic():
        return GEOGRAPHIC_SLOPE_SCALE
    return 1.0


def produce(session, product: TerrainProduct, run: Callable, *, scale=None):
    """Produce one terrain derivative, reusing an existing file.

    `run(algorithm, params)` performs the QGIS Processing call -- passed in so
    this module never imports QGIS. Returns the output path, or None on failure
    or missing input.
    """
    dem_path = session.dem_source
    if not dem_path or not os.path.exists(dem_path):
        session.tell(WARNING, "Input Error", "Input DEM file not found.")
        return None

    output = os.path.join(geometry_folder(session.project_folder), product.filename)
    setattr(session, product.field, output)

    if os.path.exists(output):
        session.say(f"{product.label} already exists. Loading existing file...")
        session.show(output, product.layer_name)
        return output

    session.say(f"Processing {product.label}...")
    if product is SLOPE:
        scale = slope_scale(session.crs) if scale is None else scale
        if scale != 1.0:
            session.say(
                f"Using gdaldem geographic scale factor for slope: {scale}")
        params = slope_params(dem_path, output, scale)
    else:
        params = aspect_params(dem_path, output)

    if not run(product.algorithm, params):
        session.say(f"{product.label} processing failed.")
        return None

    session.show(output, product.layer_name)
    session.say(f"{product.label} processing completed successfully.")
    return output


__all__ = [
    "ASPECT",
    "GEOGRAPHIC_SLOPE_SCALE",
    "SLOPE",
    "TerrainProduct",
    "aspect_params",
    "produce",
    "slope_params",
    "slope_scale",
]

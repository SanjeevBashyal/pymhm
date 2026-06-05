# -*- coding: utf-8 -*-
"""Domain discovery for meteorology forcing preparation."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from ..project_layout import geometry_folder


@dataclass(frozen=True)
class MeteoDomainBounds:
    """WGS84 bounds used to crop meteo forcing files."""

    west: float
    east: float
    south: float
    north: float
    source: str

    def as_tuple(self) -> tuple[float, float, float, float]:
        return self.west, self.east, self.south, self.north


def bounds_from_dialog(dialog) -> MeteoDomainBounds | None:
    """
    Get WGS84 crop bounds from the prepared masked DEM or selected DEM.

    Returns None when no valid domain can be discovered; callers can then keep
    the full ERA5-Land domain.
    """
    try:
        from qgis.core import (
            QgsCoordinateReferenceSystem,
            QgsCoordinateTransform,
            QgsProject,
            QgsRasterLayer,
        )
    except Exception as e:
        dialog.log_message(
            f"WARNING: QGIS domain tools unavailable; meteo data will not be cropped. {e}")
        return None

    layer = None
    source = ""

    if dialog.project_folder:
        masked_dem = Path(geometry_folder(dialog.project_folder)) / "1_dem_filled_masked.tif"
        if masked_dem.exists():
            candidate = QgsRasterLayer(str(masked_dem), "DEM_Masked")
            if candidate.isValid():
                layer = candidate
                source = str(masked_dem)

    if layer is None and hasattr(dialog, "mMapLayerComboBox_dem"):
        candidate = dialog.mMapLayerComboBox_dem.currentLayer()
        if candidate:
            layer = candidate
            source = candidate.source()

    if layer is None or not layer.isValid():
        dialog.log_message(
            "No DEM domain available for meteorology cropping; using full ERA5-Land extent.")
        return None

    layer_crs = layer.crs()
    if not layer_crs.isValid() and hasattr(dialog, "get_crs"):
        layer_crs = dialog.get_crs()

    if not layer_crs.isValid():
        dialog.log_message(
            "DEM CRS is invalid; meteo data will not be cropped.")
        return None

    try:
        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        transform = QgsCoordinateTransform(
            layer_crs, wgs84, QgsProject.instance())
        transform.setBallparkTransformsAreAppropriate(True)
        extent = transform.transformBoundingBox(layer.extent())
    except Exception as e:
        dialog.log_message(
            f"WARNING: Could not transform DEM extent to WGS84; using full ERA5-Land extent. {e}")
        return None

    if extent.isEmpty():
        dialog.log_message(
            "DEM extent is empty after transformation; using full ERA5-Land extent.")
        return None

    bounds = MeteoDomainBounds(
        west=float(extent.xMinimum()),
        east=float(extent.xMaximum()),
        south=float(extent.yMinimum()),
        north=float(extent.yMaximum()),
        source=os.path.basename(source.split("|")[0]) or source,
    )
    dialog.log_message(
        "Meteorology crop domain from "
        f"{bounds.source}: west={bounds.west:.6f}, east={bounds.east:.6f}, "
        f"south={bounds.south:.6f}, north={bounds.north:.6f}")
    return bounds

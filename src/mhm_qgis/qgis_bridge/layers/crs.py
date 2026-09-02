# -*- coding: utf-8 -*-
"""CRS questions that are answered by opening a layer."""
from __future__ import annotations

from .loader import open_layer


def crs_of(path, *, is_raster: bool = True):
    """Return the CRS of the raster or vector at `path`, or None if unreadable."""
    layer = open_layer(path, "crs_probe", is_raster=is_raster)
    if layer is None:
        return None
    crs = layer.crs()
    return crs if crs.isValid() else None


def transform_to_raster(source_layer, raster_path, log=None):
    """Return a transform from `source_layer` onto the raster's CRS.

    None when either CRS is unusable or the two already agree -- the caller
    then needs no transform at all, so None is the ordinary answer rather than
    a failure.
    """
    from qgis.core import QgsCoordinateTransform, QgsProject

    target_crs = crs_of(raster_path)
    if target_crs is None:
        return None

    source_crs = source_layer.crs()
    if not source_crs.isValid():
        return None
    if source_crs.authid() == target_crs.authid():
        return None

    transform = QgsCoordinateTransform(
        source_crs,
        target_crs,
        QgsProject.instance(),
    )
    transform.setBallparkTransformsAreAppropriate(True)
    if log is not None:
        log(f"Transforming points from {source_crs.authid()} "
            f"to raster CRS {target_crs.authid()}.")
    return transform


__all__ = ["crs_of", "transform_to_raster"]

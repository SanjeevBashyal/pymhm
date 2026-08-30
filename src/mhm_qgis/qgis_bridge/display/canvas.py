# -*- coding: utf-8 -*-
"""Put a prepared raster into the QGIS map canvas."""
from __future__ import annotations

from pathlib import Path


def raster_source(path, variable: str | None = None) -> str:
    """Return the layer source, addressing a NetCDF variable when named."""
    path = Path(path)
    if variable and path.suffix.lower() == ".nc":
        return f'NETCDF:"{path}":{variable}'
    return str(path)


def existing_layer(name: str):
    """Return the loaded layer with this name, if the project already has one."""
    try:
        from qgis.core import QgsProject
    except ImportError:
        return None
    matches = QgsProject.instance().mapLayersByName(name)
    return matches[0] if matches else None


def select_band(layer, band: int, log=None) -> bool:
    """Point a raster layer at one band, keeping the current contrast."""
    if layer is None or not band:
        return False
    try:
        from qgis.core import QgsSingleBandGrayRenderer

        renderer = layer.renderer()
        # Reusing the renderer keeps one stretch across the whole series; a
        # fresh one would re-scale per band and flatten the differences.
        if isinstance(renderer, QgsSingleBandGrayRenderer):
            renderer.setGrayBand(band)
        else:
            layer.setRenderer(QgsSingleBandGrayRenderer(layer.dataProvider(), band))
        layer.triggerRepaint()
        return True
    except Exception as error:
        if log is not None:
            log(f"WARNING: Could not select band {band}: {error}")
        return False


def show_raster(dialog, path, name, *, variable=None, band=None, is_raster=True):
    """Show one raster, reusing the loaded layer when the name already exists.

    A time scrubber calls this on every step, and `load_layer` does not
    uniquify names, so reloading would both stack duplicates in the layer tree
    and re-open the file each tick.
    """
    layer = existing_layer(name)
    if layer is None:
        layer = dialog.load_layer(raster_source(path, variable), name, is_raster)
    if layer is not None and band:
        select_band(layer, band, log=getattr(dialog, "log_message", None))
    return layer


__all__ = ["existing_layer", "raster_source", "select_band", "show_raster"]

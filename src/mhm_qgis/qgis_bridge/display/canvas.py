# -*- coding: utf-8 -*-
"""Put a prepared raster into the QGIS map canvas."""
from __future__ import annotations

import hashlib
import tempfile
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


def _materialised_path(name: str) -> Path:
    """Return a stable scratch path for one layer name.

    Reusing the same file per layer keeps the scratch directory bounded while a
    slider is scrubbed, instead of leaving one raster behind per step.
    """
    digest = hashlib.sha1(name.encode("utf-8")).hexdigest()[:12]
    folder = Path(tempfile.gettempdir()) / "mhm_qgis_display"
    folder.mkdir(parents=True, exist_ok=True)
    return folder / f"{digest}.tif"


def show_dataarray(dialog, data, name: str):
    """Show a georeferenced DataArray, reusing the loaded layer by name.

    mHM output carries no CRS, so it cannot be addressed through a
    ``NETCDF:"path":var`` URI the way prepared forcing can -- GDAL would place it
    at the origin with unit pixels. The caller attaches a CRS and the array is
    written to a scratch raster instead.
    """
    log = getattr(dialog, "log_message", None)
    if data.rio.crs is None and log is not None:
        log(f"WARNING: {name} has no CRS; the raster cannot be placed correctly.")
    target = _materialised_path(name)
    try:
        data.rio.to_raster(target)
    except Exception as error:
        if log is not None:
            log(f"ERROR: Could not prepare {name} for display: {error}")
        return None

    layer = existing_layer(name)
    if layer is None:
        return dialog.load_layer(str(target), name, True)
    # The file behind the layer was just rewritten in place.
    layer.dataProvider().reloadData()
    layer.triggerRepaint()
    return layer


__all__ = [
    "existing_layer",
    "raster_source",
    "select_band",
    "show_dataarray",
    "show_raster",
]

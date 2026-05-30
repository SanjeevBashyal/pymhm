# -*- coding: utf-8 -*-
"""Compatibility aggregate for elevation-band raster generation."""
from .elevation_band_dialog import ElevationBandDialogMixin
from .raster_masks import RasterMaskMixin
from .watershed_rasters import WatershedRasterDiscoveryMixin
from .elevation_band_rasters import ElevationBandRasterMixin


class ElevationBandGenerationMixin(
    ElevationBandDialogMixin,
    RasterMaskMixin,
    WatershedRasterDiscoveryMixin,
    ElevationBandRasterMixin,
):
    """Compatibility aggregate for elevation-band raster generation."""

    pass

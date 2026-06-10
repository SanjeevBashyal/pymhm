# -*- coding: utf-8 -*-
"""Compatibility aggregate for soil lookup and raster processing."""
from .soil_lookup import SoilLookupMixin
from .soil_raster import SoilRasterMixin


class SoilProcessingMixin(
    SoilRasterMixin,
    SoilLookupMixin,
):
    """Compatibility aggregate for soil lookup and raster processing."""

    pass

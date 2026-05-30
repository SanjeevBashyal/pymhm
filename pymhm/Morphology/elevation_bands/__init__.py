# -*- coding: utf-8 -*-
"""Compatibility aggregate for elevation-band morphology mixins."""
from .elevation_band_generation import ElevationBandGenerationMixin
from .band_details import BandDetailsMixin


class ElevationBandsMixin(
    ElevationBandGenerationMixin,
    BandDetailsMixin,
):
    """Compatibility aggregate for elevation-band morphology mixins."""

    pass

# -*- coding: utf-8 -*-
"""Compatibility aggregate for elevation-band morphology mixins."""
from .elevation_band_generation import ElevationBandGenerationMixin
from .band_details import BandDetailsMixin


class ElevationBandsMixin(
    BandDetailsMixin,
    ElevationBandGenerationMixin,
):
    """Compatibility aggregate for elevation-band morphology mixins."""

    pass

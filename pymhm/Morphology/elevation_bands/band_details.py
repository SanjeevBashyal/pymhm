# -*- coding: utf-8 -*-
"""Compatibility aggregate for elevation-band detail-table helpers."""
from .band_summary import BandSummaryMixin
from .raster_alignment import RasterAlignmentMixin
from .band_landcover_helpers import BandLandCoverHelperMixin
from .band_landcover_details import BandLandCoverDetailsMixin


class BandDetailsMixin(
    BandLandCoverDetailsMixin,
    BandSummaryMixin,
    RasterAlignmentMixin,
    BandLandCoverHelperMixin,
):
    """Compatibility aggregate for elevation-band detail-table helpers."""

    pass

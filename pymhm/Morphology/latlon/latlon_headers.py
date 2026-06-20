# -*- coding: utf-8 -*-
"""Compatibility aggregate for lat/lon header and mhm_tools processing."""
from .header_file import HeaderFileMixin
from .latlon_processing import LatLonProcessingMixin


class LatLonHeaderMixin(
    LatLonProcessingMixin,
    HeaderFileMixin,
):
    """Compatibility aggregate for lat/lon header and mhm_tools processing."""

    pass

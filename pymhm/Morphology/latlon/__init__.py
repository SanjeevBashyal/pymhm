# -*- coding: utf-8 -*-
"""Compatibility aggregate for lat/lon and ASCII export mixins."""
from .latlon_headers import LatLonHeaderMixin
from .ascii_export import AsciiExportMixin


class LatLonMixin(
    LatLonHeaderMixin,
    AsciiExportMixin,
):
    """Compatibility aggregate for lat/lon and ASCII export mixins."""

    pass

# -*- coding: utf-8 -*-
"""Compatibility aggregate for lat/lon header and NetCDF processing."""
from .header_file import HeaderFileMixin
from .latlon_processing import LatLonProcessingMixin
from .latlon_netcdf import LatLonNetcdfMixin


class LatLonHeaderMixin(
    HeaderFileMixin,
    LatLonProcessingMixin,
    LatLonNetcdfMixin,
):
    """Compatibility aggregate for lat/lon header and NetCDF processing."""

    pass

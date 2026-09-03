# -*- coding: utf-8 -*-
"""Deprecated compatibility module.

The plugin now creates latlon.nc through mhm_tools.pre.latlon.create_latlon.
This module is intentionally kept without an independent writer so stale imports
fail loudly instead of producing a second NetCDF implementation.
"""


class LatLonNetcdfMixin:
    """Deprecated compatibility stub for removed independent NetCDF writing."""

    def create_latlon_nc_file(self, latlon_folder: str) -> bool:
        raise RuntimeError(
            "Independent latlon.nc writing has been removed. "
            "Use process_lat_lon(), which calls mhm_tools.pre.latlon.create_latlon."
        )

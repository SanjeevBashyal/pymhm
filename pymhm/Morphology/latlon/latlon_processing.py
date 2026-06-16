# -*- coding: utf-8 -*-
"""Lat/lon NetCDF processing through mhm_tools."""
from __future__ import annotations

import os
import sys

from ..common import QMessageBox
from ...project_layout import data_folder, plugin_root
from ..layers.masking import MaskingMixin
from .latlon_netcdf import LatLonNetcdfMixin


class LatLonProcessingMixin(MaskingMixin, LatLonNetcdfMixin):
    """Lat/lon grid metadata processing for mHM."""

    def process_lat_lon(self) -> bool:
        """Create data/latlon.nc from L0, L1, L11, and L2 grid headers."""
        self.log_message("\n--- Processing Lat/Lon NetCDF ---")

        if not self.check_prerequisites():
            return False

        try:
            headers = self.dialog.grid_level_headers()
        except Exception as e:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                str(e),
            )
            self.log_message(f"ERROR: Cannot prepare latlon headers: {e}")
            return False

        crs = self.dialog.get_crs()
        crs_string = None
        if crs is not None and crs.isValid():
            crs_string = crs.authid() or (
                f"EPSG:{crs.postgisSrid()}" if crs.postgisSrid() else None
            )

        output_folder = data_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, "latlon.nc")

        try:
            self._ensure_local_mhm_tools_importable()
            from mhm_tools.pre.latlon import create_latlon

            level0 = headers["L0"]
            level1 = headers["L1"]["cellsize"]
            level11 = headers["L11"]["cellsize"]
            level2 = headers["L2"]["cellsize"]
            create_latlon(
                out_file=output_file,
                level0=level0,
                level1=level1,
                level11=level11,
                level2=level2,
                crs=crs_string,
            )
        except ImportError as e:
            self.log_message(f"ERROR: Required library not available: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"Required library not available: {e}\n"
                "Please install xarray, netCDF4, numpy, and pyproj in the QGIS Python environment.",
            )
            return False
        except Exception as e:
            import traceback

            details = traceback.format_exc()
            self.log_message(f"ERROR creating latlon.nc file: {e}\n{details}")
            QMessageBox.critical(
                self.dialog,
                "Lat/Lon Error",
                f"Error creating latlon.nc:\n{e}",
            )
            return False

        self.mark_output_prepared(
            output_file,
            name="latlon.nc",
            loaded=False,
        )
        self.log_message(f"latlon.nc created successfully: {output_file}")
        return True

    def _ensure_local_mhm_tools_importable(self) -> None:
        """Put the plugin root on sys.path so bundled mhm_tools imports work."""
        root = plugin_root()
        if root not in sys.path:
            sys.path.insert(0, root)

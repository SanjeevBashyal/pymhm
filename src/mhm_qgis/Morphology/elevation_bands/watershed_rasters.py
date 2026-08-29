# -*- coding: utf-8 -*-
"""Per-subcatchment watershed raster discovery."""
from ..common import os


class WatershedRasterDiscoveryMixin:
    """Per-subcatchment watershed raster discovery."""

    def _collect_watershed_rasters(self, watershed_output_folder):
        """Find per-subcatchment watershed rasters."""
        if not os.path.exists(watershed_output_folder):
            return []

        raster_paths = []
        for filename in os.listdir(watershed_output_folder):
            lower_name = filename.lower()
            if not lower_name.endswith(".tif"):
                continue
            if not lower_name.startswith("4_watershed_"):
                continue
            if lower_name.endswith("_raw.tif"):
                continue
            raster_paths.append(os.path.join(watershed_output_folder, filename))

        return sorted(raster_paths)

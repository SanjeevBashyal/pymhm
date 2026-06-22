# -*- coding: utf-8 -*-
"""Terrain derivative outputs: aspect and slope."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
)
from ..core.dem_inputs import DemInputMixin
from ...mhm_tools_to_integrate.setup_creation.terrain import aspect_params, slope_params


class TerrainAnalysisMixin(DemInputMixin):
    """Terrain derivative outputs: aspect and slope."""

    def process_aspect(self) -> None:
        """Process Aspect from input DEM"""
        # Ensure DEM is checked and reprojected if needed
        if not self.dem_layer:
            input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            if not input_dem_layer:
                QMessageBox.warning(self.dialog, "Input Error",
                                    "Please select an input DEM layer.")
                return
            self.check_and_reproject_dem(input_dem_layer)

        # Use reprojected DEM layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()
        input_dem_path = dem_layer.source()

        if not os.path.exists(input_dem_path):
            QMessageBox.warning(self.dialog, "Input Error",
                                "Input DEM file not found.")
            return

        # Check if aspect already exists
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.aspect_path = os.path.join(
            geometry_folder, "1_dem_aspect.tif")
        if self.aspect_path and os.path.exists(self.aspect_path):
            self.log_message("Aspect already exists. Loading existing file...")
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            return

        self.log_message("Processing Aspect...")

        params_aspect = aspect_params(input_dem_path, self.aspect_path)

        result = self.run_processing_algorithm("gdal:aspect", params_aspect)
        if result:
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            self.log_message("Aspect processing completed successfully.")
        else:
            self.log_message("Aspect processing failed.")

    def process_slope(self) -> None:
        """Process Slope from input DEM"""
        # Ensure DEM is checked and reprojected if needed
        if not self.dem_layer:
            input_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            if not input_dem_layer:
                QMessageBox.warning(self.dialog, "Input Error",
                                    "Please select an input DEM layer.")
                return
            self.check_and_reproject_dem(input_dem_layer)

        # Use reprojected DEM layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()
        input_dem_path = dem_layer.source()

        if not os.path.exists(input_dem_path):
            QMessageBox.warning(self.dialog, "Input Error",
                                "Input DEM file not found.")
            return

        # Check if slope already exists
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.slope_path = os.path.join(
            geometry_folder, "1_dem_slope.tif")
        if self.slope_path and os.path.exists(self.slope_path):
            self.log_message("Slope already exists. Loading existing file...")
            self.load_layer(self.slope_path, "1_DEM_Slope")
            return

        self.log_message("Processing Slope...")

        slope_scale = self._slope_scale_for_dem(dem_layer)
        if slope_scale != 1.0:
            self.log_message(
                f"Using gdaldem geographic scale factor for slope: {slope_scale}")

        params_slope = slope_params(input_dem_path, self.slope_path, slope_scale)

        result = self.run_processing_algorithm("gdal:slope", params_slope)
        if result:
            self.load_layer(self.slope_path, "1_DEM_Slope")
            self.log_message("Slope processing completed successfully.")
        else:
            self.log_message("Slope processing failed.")

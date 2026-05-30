# -*- coding: utf-8 -*-
"""Terrain derivative outputs: aspect and slope."""
from ..common import (
    os,
    QMessageBox,
    processing,
)


class TerrainAnalysisMixin:
    """Terrain derivative outputs: aspect and slope."""

    def process_aspect(self):
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
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.aspect_path = os.path.join(
            geometry_folder, "1_dem_aspect.tif")
        if self.aspect_path and os.path.exists(self.aspect_path):
            self.log_message("Aspect already exists. Loading existing file...")
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            return

        self.log_message("Processing Aspect...")

        params_aspect = {
            'INPUT': input_dem_path,
            'BAND': 1,
            'TRIG_ANGLE': False,
            'ZERO_FLAT': False,
            'COMPUTE_EDGES': False,
            'ZEVENBERGEN': False,
            'OPTIONS': None,
            'EXTRA': '',
            'OUTPUT': self.aspect_path
        }

        result = self.run_processing_algorithm("gdal:aspect", params_aspect)
        if result:
            self.load_layer(self.aspect_path, "1_DEM_Aspect")
            self.log_message("Aspect processing completed successfully.")
        else:
            self.log_message("Aspect processing failed.")

    def process_slope(self):
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
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
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

        params_slope = {
            'INPUT': input_dem_path,
            'BAND': 1,
            'SCALE': slope_scale,
            'AS_PERCENT': True,
            'COMPUTE_EDGES': False,
            'ZEVENBERGEN': False,
            'OPTIONS': None,
            'EXTRA': '',
            'OUTPUT': self.slope_path
        }

        result = self.run_processing_algorithm("gdal:slope", params_slope)
        if result:
            self.load_layer(self.slope_path, "1_DEM_Slope")
            self.log_message("Slope processing completed successfully.")
        else:
            self.log_message("Slope processing failed.")

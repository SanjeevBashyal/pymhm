# -*- coding: utf-8 -*-
"""
Utility functions and base classes for pymhm dialog processing
"""
import os
from qgis.PyQt.QtWidgets import QMessageBox
from qgis.core import (
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsApplication
)
import processing


class DialogUtils:
    """
    Utility class providing common utility methods for dialog processing.
    Methods expect the dialog instance to have:
    - self.LogText (QTextBrowser)
    - self.project_folder (str)
    - self.geometry_folder (str)
    - self.mMapLayerComboBox_dem (QgsMapLayerComboBox)
    - self.mMapLayerComboBox_pour_points (QgsMapLayerComboBox)
    """
    
    def log_message(self, message):
        """Appends a message to the log text browser."""
        self.LogText.append(message)
        QgsApplication.processEvents()

    def check_prerequisites(self, needs_pour_points=False):
        """Check if project folder and necessary layers are set."""
        if not self.project_folder:
            QMessageBox.critical(
                self, "Missing Input", "Please select a project folder before proceeding.")
            return False
        if not self.mMapLayerComboBox_dem.currentLayer():
            QMessageBox.critical(self, "Missing Input",
                                 "Please select a DEM Raster Layer.")
            return False
        if needs_pour_points and not self.mMapLayerComboBox_pour_points.currentLayer():
            QMessageBox.critical(self, "Missing Input",
                                 "Please select a Pour Points Layer.")
            return False
        return True

    def run_processing_algorithm(self, name, params):
        """A wrapper to run a processing algorithm and handle errors."""
        self.log_message(f"Running algorithm: {name}...")
        try:
            result = processing.run(name, params)
            self.log_message(f"Algorithm '{name}' finished successfully.")
            return result
        except Exception as e:
            self.log_message(
                f"ERROR: Algorithm '{name}' failed. Details: {str(e)}")
            QMessageBox.critical(
                self, "Processing Error", f"Algorithm '{name}' failed.\nCheck the log for details.")
            return None

    def load_layer(self, path, name, is_raster=True):
        """Loads a layer into the QGIS project."""
        if not path or not os.path.exists(path):
            self.log_message(f"ERROR: Output file not found at {path}")
            return

        layer = QgsRasterLayer(
            path, name) if is_raster else QgsVectorLayer(path, name, "ogr")

        if not layer.isValid():
            self.log_message(f"ERROR: Failed to load layer: {name}")
            return

        QgsProject.instance().addMapLayer(layer)
        self.log_message(f"Layer '{name}' added to project.")

    def get_dem_extent_and_resolution(self):
        """Get DEM extent and pixel resolution for clipping and rasterization"""
        dem_layer = self.mMapLayerComboBox_dem.currentLayer()
        if not dem_layer:
            return None, None, None
        
        # Get extent
        extent = dem_layer.extent()
        extent_str = f"{extent.xMinimum()},{extent.xMaximum()},{extent.yMinimum()},{extent.yMaximum()}"
        
        # Get pixel size
        raster_extent = dem_layer.extent()
        width = dem_layer.width()
        height = dem_layer.height()
        pixel_size_x = (raster_extent.xMaximum() - raster_extent.xMinimum()) / width
        pixel_size_y = (raster_extent.yMaximum() - raster_extent.yMinimum()) / height
        
        return extent_str, pixel_size_x, pixel_size_y


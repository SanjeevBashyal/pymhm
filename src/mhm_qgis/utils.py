# -*- coding: utf-8 -*-
"""
Utility functions and base classes for mhm_qgis dialog processing
"""
import sys
from qgis.PyQt.QtWidgets import QMessageBox
from qgis.core import QgsApplication


class DialogUtils:
    """
    Utility class providing common utility methods for dialog processing.
    Methods expect the dialog instance to have:
    - self.LogText (QTextBrowser)
    - self.project_folder (str)
    - self.geometry_folder (str)
    - self.input_combo("dem") (QgsMapLayerComboBox)
    - self.input_combo("pour_points") (QgsMapLayerComboBox)
    """
    
    def log_message(self, message):
        """Appends a message to the log text browser."""
        try:
            console = getattr(sys, "__stdout__", None)
            if console is not None:
                print(message, file=console)
        except Exception:
            pass
        self.LogText.append(message)
        QgsApplication.processEvents()

    def check_prerequisites(self, needs_pour_points=False):
        """Check if project folder and necessary layers are set."""
        if not self.project_folder:
            QMessageBox.critical(
                self, "Missing Input", "Please select a project folder before proceeding.")
            return False
        if not self.input_combo("dem").currentLayer():
            QMessageBox.critical(self, "Missing Input",
                                 "Please select a DEM Raster Layer.")
            return False
        if needs_pour_points and not self.input_combo("pour_points").currentLayer():
            QMessageBox.critical(self, "Missing Input",
                                 "Please select a Pour Points Layer.")
            return False
        return True

    def run_processing_algorithm(self, name, params):
        """Run a processing algorithm, reporting failures to the user."""
        from .qgis_bridge import processing as qgis_processing
        from .qt.reporting import dialog_reporter

        return qgis_processing.run(
            name, params,
            log=self.log_message,
            report=dialog_reporter(self, log=None),
        )

    def get_dem_extent_and_resolution(self):
        """Get DEM extent and pixel resolution for clipping and rasterization"""
        # Check if morphology processor has reprojected DEM layer
        if hasattr(self, 'morphology_processor') and self.morphology_processor.dem_layer:
            dem_layer = self.morphology_processor.dem_layer
        else:
            dem_layer = self.input_combo("dem").currentLayer()
        
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

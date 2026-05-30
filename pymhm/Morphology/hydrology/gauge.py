# -*- coding: utf-8 -*-
"""Gauge-position rasterization from snapped pour points."""
from ..common import (
    os,
    QMessageBox,
    QgsRasterLayer,
    processing,
)


class GaugePositionMixin:
    """Gauge-position rasterization from snapped pour points."""

    def process_gauge_position(self):
        """Process Gauge Position from pour points"""
        # Check prerequisites
        if not self._ensure_filled_dem():
            return

        if not self.snapped_points_path or not os.path.exists(self.snapped_points_path):
            if not self._ensure_snapped_points():
                return

        # Check if gauge position already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.gauge_position_path = os.path.join(
            geometry_folder, "2_gauge_position.tif")
        if self.gauge_position_path and os.path.exists(self.gauge_position_path):
            self.log_message(
                "Gauge Position already exists. Loading existing file...")
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            return

        self.log_message("Processing Gauge Position...")

        # Get cell size from filled DEM
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM to get cell size.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM to get cell size.")
            return

        # Calculate cell size (width and height) and get extent
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        cell_size_width = (raster_extent.xMaximum() -
                           raster_extent.xMinimum()) / width
        cell_size_height = (raster_extent.yMaximum() -
                            raster_extent.yMinimum()) / height

        # Format extent as string: "xmin,xmax,ymin,ymax"
        extent_str = f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},{raster_extent.yMinimum()},{raster_extent.yMaximum()}"

        self.log_message(
            f"Cell size from filled DEM - Width: {cell_size_width:.2f}m, Height: {cell_size_height:.2f}m")
        self.log_message(f"Using extent from filled DEM: {extent_str}")

        params_gauge = {
            'INPUT': self.snapped_points_path,
            'FIELD': 'Name',
            'BURN': 0,
            'USE_Z': False,
            'UNITS': 1,
            'WIDTH': cell_size_width,
            'HEIGHT': cell_size_height,
            'EXTENT': extent_str,
            'NODATA': -9999,
            'OPTIONS': None,
            'DATA_TYPE': 5,
            'INIT': None,
            'INVERT': False,
            'EXTRA': '',
            'OUTPUT': self.gauge_position_path
        }

        result = self.run_processing_algorithm("gdal:rasterize", params_gauge)
        if result:
            self.load_layer(self.gauge_position_path, "2_Gauge_Position")
            self.log_message(
                "Gauge Position processing completed successfully.")
        else:
            self.log_message("Gauge Position processing failed.")

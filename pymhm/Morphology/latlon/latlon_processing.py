# -*- coding: utf-8 -*-
"""Lat/lon grid metadata processing for mHM."""
from ..common import (
    os,
    QMessageBox,
    QgsRasterLayer,
    geometry_folder,
    processing,
)
from ...project_layout import data_folder


class LatLonProcessingMixin:
    """Lat/lon grid metadata processing for mHM."""

    def process_lat_lon(self):
        """
        Process lat/lon information by reading dem_masked_layer and creating header files
        for L0, L1, L11, and L2 levels. This prepares data for latlon.nc file creation.
        """
        self.log_message("\n--- Processing Lat/Lon Headers ---")
        
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        geom_folder = geometry_folder(self.dialog.project_folder)
        dem_masked_path = os.path.join(geom_folder, "1_dem_filled_masked.tif")
        
        if not os.path.exists(dem_masked_path):
            self.log_message("Masked DEM is missing. Running Mask All Layers first...")
            self.without_layer_loading(self.mask_all_layers)
            if not os.path.exists(dem_masked_path):
                return
        
        # Read L0, L1, L2 values and units from UI
        try:
            l0_value = float(self.dialog.lineEdit_L0.text())
            l0_unit = self.dialog.comboBox_L0.currentText()
            l1_value = float(self.dialog.lineEdit_L1.text())
            l1_unit = self.dialog.comboBox_L1.currentText()
            l2_value = float(self.dialog.lineEdit_L2.text())
            l2_unit = self.dialog.comboBox_L2.currentText()
        except (ValueError, AttributeError) as e:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Please enter valid numeric values for L0, L1, and L2.\nError: {e}")
            return
        
        if not l0_unit or not l1_unit or not l2_unit:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select units (m or °) for L0, L1, and L2.")
            return
        
        # Read the masked DEM layer
        dem_masked_layer = QgsRasterLayer(dem_masked_path, "DEM_Masked")
        if not dem_masked_layer.isValid():
            self.log_message("ERROR: Cannot read masked DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read masked DEM layer.")
            return
        
        # Get extent and dimensions from masked DEM
        raster_extent = dem_masked_layer.extent()
        width = dem_masked_layer.width()  # Number of columns
        height = dem_masked_layer.height()  # Number of rows
        
        # Calculate cell size (assuming square cells)
        cell_size_x = (raster_extent.xMaximum() - raster_extent.xMinimum()) / width
        cell_size_y = (raster_extent.yMaximum() - raster_extent.yMinimum()) / height
        cell_size = (cell_size_x + cell_size_y) / 2.0  # Average cell size
        
        # Get corner coordinates
        xllcorner = raster_extent.xMinimum()
        yllcorner = raster_extent.yMinimum()
        
        # L0: Use values from masked DEM
        self.L0 = {
            'ncols': width,
            'nrows': height,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': cell_size,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L0 information extracted from masked DEM:")
        self.log_message(f"  ncols: {self.L0['ncols']}, nrows: {self.L0['nrows']}")
        self.log_message(f"  xllcorner: {self.L0['xllcorner']}, yllcorner: {self.L0['yllcorner']}")
        self.log_message(f"  cellsize: {self.L0['cellsize']:.2f}")
        
        # Convert L1 cell size to meters if needed
        l1_cellsize = l1_value
        if l1_unit == "°":
            # Convert degrees to meters (approximate: 1 degree ≈ 111,000 meters)
            center_lat = (raster_extent.yMinimum() + raster_extent.yMaximum()) / 2.0
            l1_cellsize = l1_value * 111000
            self.log_message(f"L1 cell size converted from {l1_value}° to {l1_cellsize:.2f} m")
        
        # L1: Calculate new ncols and nrows based on L1 cell size
        # Keep same xllcorner and yllcorner
        x_range = raster_extent.xMaximum() - raster_extent.xMinimum()
        y_range = raster_extent.yMaximum() - raster_extent.yMinimum()
        
        l1_ncols = int(round(x_range / l1_cellsize))
        l1_nrows = int(round(y_range / l1_cellsize))
        
        self.L1 = {
            'ncols': l1_ncols,
            'nrows': l1_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l1_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L1 information calculated:")
        self.log_message(f"  ncols: {self.L1['ncols']}, nrows: {self.L1['nrows']}")
        self.log_message(f"  xllcorner: {self.L1['xllcorner']}, yllcorner: {self.L1['yllcorner']}")
        self.log_message(f"  cellsize: {self.L1['cellsize']:.2f}")
        
        # L11: Exactly same as L1
        self.L11 = {
            'ncols': l1_ncols,
            'nrows': l1_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l1_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L11 information (same as L1):")
        self.log_message(f"  ncols: {self.L11['ncols']}, nrows: {self.L11['nrows']}")
        self.log_message(f"  cellsize: {self.L11['cellsize']:.2f}")
        
        # Convert L2 cell size to meters if needed
        l2_cellsize = l2_value
        if l2_unit == "°":
            # Convert degrees to meters
            l2_cellsize = l2_value * 111000
            self.log_message(f"L2 cell size converted from {l2_value}° to {l2_cellsize:.2f} m")
        
        # L2: Calculate new ncols and nrows based on L2 cell size
        l2_ncols = int(round(x_range / l2_cellsize))
        l2_nrows = int(round(y_range / l2_cellsize))
        
        self.L2 = {
            'ncols': l2_ncols,
            'nrows': l2_nrows,
            'xllcorner': xllcorner,
            'yllcorner': yllcorner,
            'cellsize': l2_cellsize,
            'NODATA_value': -9999
        }
        
        self.log_message(f"L2 information calculated:")
        self.log_message(f"  ncols: {self.L2['ncols']}, nrows: {self.L2['nrows']}")
        self.log_message(f"  xllcorner: {self.L2['xllcorner']}, yllcorner: {self.L2['yllcorner']}")
        self.log_message(f"  cellsize: {self.L2['cellsize']:.2f}")
        
        latlon_folder = data_folder(self.dialog.project_folder)
        os.makedirs(latlon_folder, exist_ok=True)
        
        # Create NetCDF file
        self.log_message("\nCreating latlon.nc file...")
        self.create_latlon_nc_file(latlon_folder)
        
        self.log_message("Lat/Lon header processing completed successfully.")
        QMessageBox.information(
            self.dialog, "Success",
            "Lat/Lon headers and NetCDF file created successfully in data folder.")

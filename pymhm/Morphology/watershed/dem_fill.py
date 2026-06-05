# -*- coding: utf-8 -*-
"""DEM depression filling workflow."""
from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
)


class DemFillMixin:
    """DEM depression filling workflow."""

    def fill_dem(self):
        """Step 1: Fill sinks in the DEM using pyflwdir, as in the workshop notebook."""
        self.log_message("\n--- Starting Geometry Step 1: Fill DEM ---")
        if not self.check_prerequisites():
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.filled_dem_path = os.path.join(
            geometry_folder, "1_dem_filled.tif")

        # Check if filled DEM already exists
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            self.log_message(
                "Filled DEM intermediate already exists. Using existing file.")
            return

        original_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        self.log_message(f"Input DEM: {original_dem_layer.name()}")

        # Check and reproject DEM if needed (stores in self.dem_layer)
        dem_layer = self.check_and_reproject_dem(original_dem_layer)

        # Ensure self.dem_layer is set
        if not self.dem_layer:
            self.dem_layer = dem_layer

        deps = self._get_python_morphology_deps()
        if not deps:
            self.filled_dem_path = None
            return

        dem_source = self.dem_layer.source()
        reference = self._read_raster_array(dem_source, as_float=True)
        if not reference:
            self.filled_dem_path = None
            return

        dem, invalid_mask, pyflwdir_nodata = self._normalise_dem_array(
            reference["array"], reference["nodata"])
        if dem is None:
            self.filled_dem_path = None
            return

        self.log_message("Filling DEM depressions with pyflwdir.dem.fill_depressions...")
        try:
            filled_dem, filled_mask = deps["pfd"].dem.fill_depressions(
                dem, nodata=pyflwdir_nodata)
        except Exception as e:
            self.log_message(f"ERROR: pyflwdir DEM filling failed: {e}")
            QMessageBox.critical(
                self.dialog, "Processing Error",
                f"pyflwdir DEM filling failed.\n{e}")
            self.filled_dem_path = None
            return

        filled_dem = filled_dem.astype(deps["np"].float32)
        filled_dem[invalid_mask] = pyflwdir_nodata

        if self._write_raster_array(
                self.filled_dem_path,
                filled_dem,
                reference,
                nodata=pyflwdir_nodata,
                gdal_type=deps["gdal"].GDT_Float32):
            filled_cells = int(deps["np"].count_nonzero(
                filled_mask)) if filled_mask is not None else 0
            self.log_message(
                f"DEM filled successfully with pyflwdir ({filled_cells} adjusted cells).")
        else:
            self.log_message("ERROR: Failed to write filled DEM raster.")
            self.filled_dem_path = None

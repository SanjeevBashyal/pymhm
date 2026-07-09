# -*- coding: utf-8 -*-
"""DEM depression filling workflow."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsRasterLayer,
)
from ..core.dem_inputs import DemInputMixin
from ..core.raster_io import RasterIOMixin


class DemFillMixin(DemInputMixin, RasterIOMixin):
    """DEM depression filling workflow."""

    def fill_dem(self) -> None:
        """Step 1: Fill sinks in the DEM using pyflwdir, as in the workshop notebook."""
        self.log_message("\n--- Starting Geometry Step 1: Fill DEM ---")
        if not self.check_prerequisites():
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        self.filled_dem_path = os.path.join(
            geometry_folder, "1_dem_filled.tif")

        original_dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        self.log_message(f"Input DEM: {original_dem_layer.name()}")

        # Check and reproject DEM if needed (stores in self.dem_layer)
        dem_layer = self.check_and_reproject_dem(original_dem_layer)

        # Ensure self.dem_layer is set
        if not self.dem_layer:
            self.dem_layer = dem_layer

        # Check if filled DEM already exists and still matches the current DEM grid.
        if self.filled_dem_path and os.path.exists(self.filled_dem_path):
            if self._existing_filled_dem_matches_current_dem(
                    self.filled_dem_path, self.dem_layer):
                self.log_message(
                    "Filled DEM intermediate already exists. Using existing file.")
                self.load_layer(self.filled_dem_path, "1_DEM_Filled")
                return
            self.log_message(
                "Existing filled DEM does not match the current DEM CRS/grid. "
                "Recreating dependent geometry outputs.")
            self._remove_stale_flow_outputs()

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
            self.load_layer(self.filled_dem_path, "1_DEM_Filled")
            self.log_message(
                f"DEM filled successfully with pyflwdir ({filled_cells} adjusted cells).")
        else:
            self.log_message("ERROR: Failed to write filled DEM raster.")
            self.filled_dem_path = None

    def _existing_filled_dem_matches_current_dem(
            self,
            filled_dem_path: str,
            dem_layer: object) -> bool:
        """Return True when an existing filled DEM matches the active DEM grid."""
        if not dem_layer or not os.path.exists(filled_dem_path):
            return False

        filled_layer = QgsRasterLayer(filled_dem_path, "Existing_Filled_DEM")
        if not filled_layer.isValid():
            return False

        if filled_layer.width() != dem_layer.width():
            return False
        if filled_layer.height() != dem_layer.height():
            return False

        filled_crs = filled_layer.crs()
        dem_crs = dem_layer.crs()
        if filled_crs.isValid() != dem_crs.isValid():
            return False
        if filled_crs.isValid() and filled_crs.authid() != dem_crs.authid():
            return False

        tolerance = 1e-9
        for filled_value, dem_value in (
                (filled_layer.rasterUnitsPerPixelX(), dem_layer.rasterUnitsPerPixelX()),
                (filled_layer.rasterUnitsPerPixelY(), dem_layer.rasterUnitsPerPixelY()),
                (filled_layer.extent().xMinimum(), dem_layer.extent().xMinimum()),
                (filled_layer.extent().xMaximum(), dem_layer.extent().xMaximum()),
                (filled_layer.extent().yMinimum(), dem_layer.extent().yMinimum()),
                (filled_layer.extent().yMaximum(), dem_layer.extent().yMaximum())):
            scale = max(abs(float(filled_value)), abs(float(dem_value)), 1.0)
            if abs(float(filled_value) - float(dem_value)) > tolerance * scale:
                return False

        return True

    def _remove_stale_flow_outputs(self) -> None:
        """Remove geometry outputs that depend on the filled DEM grid."""
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        raster_names = (
            "1_dem_filled.tif",
            "1_dem_filled.tif.aux.xml",
            "2_flow_accumulation.tif",
            "2_flow_accumulation.tif.aux.xml",
            "2_flow_accumulation_area.tif",
            "2_flow_accumulation_area.tif.aux.xml",
            "2_flow_direction.tif",
            "2_flow_direction.tif.aux.xml",
            "2_gauge_position.tif",
            "2_gauge_position.tif.aux.xml",
        )
        for filename in raster_names:
            path = os.path.join(geometry_folder, filename)
            if os.path.exists(path):
                try:
                    os.remove(path)
                except Exception as e:
                    self.log_message(f"WARNING: Could not remove stale file {path}: {e}")

        for vector_name in (
                "2_channel_network.shp",
                "2_pour_points_snapped.shp"):
            vector_path = os.path.join(geometry_folder, vector_name)
            if hasattr(self, "_remove_vector_dataset"):
                self._remove_vector_dataset(vector_path)

        watersheds_folder = os.path.join(geometry_folder, "Watersheds")
        if os.path.isdir(watersheds_folder):
            for root, _dirs, files in os.walk(watersheds_folder):
                for filename in files:
                    path = os.path.join(root, filename)
                    try:
                        os.remove(path)
                    except Exception as e:
                        self.log_message(f"WARNING: Could not remove stale file {path}: {e}")

        self.flow_direction_path = None
        self.flow_accumulation_path = None
        self.flow_accumulation_area_path = None
        self.channel_network_vector_path = None
        self.snapped_points_path = None
        self.gauge_position_path = None
        self.watershed_raster_path = None
        self.watershed_vector_path = None
        self.merged_watershed_path = None

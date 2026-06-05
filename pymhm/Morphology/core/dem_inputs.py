# -*- coding: utf-8 -*-
"""DEM validation, reprojection, and ASCII preparation."""
from ..common import (
    os,
    QMessageBox,
    QgsRasterLayer,
    geometry_folder,
    morph_folder,
    processing,
)


class DemInputMixin:
    """DEM validation, reprojection, and ASCII preparation."""

    def check_and_reproject_dem(self, dem_layer):
        """
        Check if DEM CRS matches input CRS, and reproject if needed.
        Stores the reprojected layer in self.dem_layer attribute.

        Args:
            dem_layer: The DEM raster layer

        Returns:
            The DEM layer (original or reprojected) to use for processing
        """
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            self.log_message(
                "WARNING: Input CRS is not valid. Using DEM CRS as-is.")
            self.dem_layer = dem_layer
            return dem_layer

        dem_crs = dem_layer.crs()
        if not dem_crs.isValid():
            self.log_message("WARNING: DEM CRS is not valid. Using DEM as-is.")
            self.dem_layer = dem_layer
            return dem_layer

        # Check if CRS matches
        if dem_crs.authid() == input_crs.authid():
            self.log_message(
                f"DEM CRS ({dem_crs.authid()}) matches input CRS. No reprojection needed.")
            self.dem_layer = dem_layer
            return dem_layer

        # Need to reproject - check if reprojected file already exists
        reprojected_dem_path = os.path.join(
            geometry_folder(self.dialog.project_folder),
            "0_dem_reprojected.tif")

        if os.path.exists(reprojected_dem_path):
            self.log_message(f"Found existing reprojected DEM. Loading it...")
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message(
                    "WARNING: Existing reprojected DEM is not valid. Reprojecting again...")

        # Reproject DEM
        self.log_message(
            f"DEM CRS ({dem_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")

        params_warp = {
            'INPUT': dem_layer.source(),
            'SOURCE_CRS': None,  # Auto-detect from input
            'TARGET_CRS': input_crs,
            'RESAMPLING': 0,  # Nearest neighbor
            'NODATA': None,
            'TARGET_RESOLUTION': None,
            'OPTIONS': None,
            'DATA_TYPE': 0,  # Use input data type
            'TARGET_EXTENT': None,
            'TARGET_EXTENT_CRS': None,
            'MULTITHREADING': False,
            'EXTRA': '',
            'OUTPUT': reprojected_dem_path
        }

        result = self.run_processing_algorithm(
            "gdal:warpreproject", params_warp)
        if result and os.path.exists(reprojected_dem_path):
            reprojected_layer = QgsRasterLayer(
                reprojected_dem_path, "1_DEM_Reprojected")
            if reprojected_layer.isValid():
                self.log_message(
                    f"DEM reprojected successfully to {input_crs.authid()}")
                # Do not add reprojected layer to QGIS project - it's only used for processing
                self.dem_layer = reprojected_layer
                return reprojected_layer
            else:
                self.log_message(
                    "ERROR: Reprojected DEM layer is not valid. Using original DEM.")
                self.dem_layer = dem_layer
                return dem_layer
        else:
            self.log_message(
                "ERROR: DEM reprojection failed. Using original DEM.")
            self.dem_layer = dem_layer
            return dem_layer

    def convert_dem_to_asc(self):
        """Convert DEM from TIF to ASC format. First reprojects to input CRS, then converts to ASC."""
        self.log_message("\n--- Converting DEM to ASC format ---")
        if not self.check_prerequisites():
            return

        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        dem_asc_path = os.path.join(morph_output_folder, "dem.asc")
        
        # Geometry folder for temporary files
        geom_folder = geometry_folder(self.dialog.project_folder)

        # Check if ASC file already exists
        if os.path.exists(dem_asc_path):
            self.log_message(
                "DEM ASC file already exists. Loading existing file...")
            self.load_layer(dem_asc_path, "DEM_ASC", is_raster=True)
            return

        # Get DEM layer from combo box
        dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not dem_layer:
            QMessageBox.critical(
                self.dialog, "Missing Input",
                "Please select a DEM Raster Layer.")
            return

        # Get input CRS from CRS widget
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.critical(
                self.dialog, "Missing Input",
                "Please select a valid Coordinate Reference System.")
            return

        self.log_message(f"Input DEM: {dem_layer.name()}")

        # Get DEM CRS
        dem_crs = dem_layer.crs()
        if not dem_crs.isValid():
            self.log_message(
                "WARNING: DEM CRS is not valid. Proceeding with conversion...")

        self.log_message(
            f"DEM CRS: {dem_crs.authid() if dem_crs.isValid() else 'Unknown'}")
        self.log_message(f"Target CRS: {input_crs.authid()}")

        # Get DEM source path
        dem_source = dem_layer.source()
        if not dem_source or not os.path.exists(dem_source):
            self.log_message(f"ERROR: DEM source file not found: {dem_source}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"DEM source file not found: {dem_source}")
            return

        # Step 1: Reproject DEM to input CRS (if needed)
        reprojected_dem_path = os.path.join(
            geom_folder, "dem_reprojected.tif")
        dem_source_for_conversion = dem_source

        # Check if reprojection is needed
        needs_reprojection = False
        if dem_crs.isValid() and input_crs.isValid():
            if dem_crs.authid() != input_crs.authid():
                needs_reprojection = True

        if needs_reprojection:
            self.log_message(
                f"Reprojecting DEM from {dem_crs.authid()} to {input_crs.authid()}...")

            # Check if reprojected file already exists
            if os.path.exists(reprojected_dem_path):
                self.log_message("Found existing reprojected DEM. Using it...")
                dem_source_for_conversion = reprojected_dem_path
            else:
                # Reproject DEM
                params_warp = {
                    'INPUT': dem_source,
                    'SOURCE_CRS': None,  # Auto-detect from input
                    'TARGET_CRS': input_crs,
                    'RESAMPLING': 0,  # Nearest neighbor
                    'NODATA': None,
                    'TARGET_RESOLUTION': None,
                    'OPTIONS': None,
                    'DATA_TYPE': 0,  # Use input data type
                    'TARGET_EXTENT': None,
                    'TARGET_EXTENT_CRS': None,
                    'MULTITHREADING': False,
                    'EXTRA': '',
                    'OUTPUT': reprojected_dem_path
                }

                result_warp = self.run_processing_algorithm(
                    "gdal:warpreproject", params_warp)
                if result_warp and os.path.exists(reprojected_dem_path):
                    self.log_message("DEM reprojected successfully!")
                    dem_source_for_conversion = reprojected_dem_path
                else:
                    self.log_message(
                        "ERROR: DEM reprojection failed. Using original DEM.")
                    QMessageBox.warning(
                        self.dialog, "Warning",
                        "DEM reprojection failed. Converting original DEM to ASC.")
        else:
            self.log_message(
                "DEM CRS matches target CRS. No reprojection needed.")

        # Step 2: Convert to ASC format
        self.log_message(f"Converting DEM to ASC format...")
        self.log_message(f"Output file: {dem_asc_path}")

        # Parameters for GDAL translate
        params_translate = {
            'INPUT': dem_source_for_conversion,
            'TARGET_CRS': None,  # Already in target CRS from reprojection
            'NODATA': -9999,
            'COPY_SUBDATASETS': False,
            'OPTIONS': None,
            'EXTRA': '',
            'DATA_TYPE': 0,  # Use input data type
            'OUTPUT': dem_asc_path
        }

        result = self.run_processing_algorithm(
            "gdal:translate", params_translate)
        if result and os.path.exists(dem_asc_path):
            self.load_layer(dem_asc_path, "DEM_ASC", is_raster=True)
            self.log_message("DEM converted to ASC format successfully!")
        else:
            self.log_message("ERROR: Failed to convert DEM to ASC format.")

    def _slope_scale_for_dem(self, dem_layer):
        """Match the notebook's gdaldem scale choice for geographic DEMs."""
        crs = dem_layer.crs()
        if crs.isValid() and crs.isGeographic():
            return 111200.0
        return 1.0

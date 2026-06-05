# -*- coding: utf-8 -*-
"""Generic layer reprojection, clipping, and rasterization helpers."""
from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsVectorLayer,
    QgsRasterLayer,
    QgsCoordinateReferenceSystem,
    processing,
)


class LayerPreparationMixin:
    """Generic layer reprojection, clipping, and rasterization helpers."""

    def reproject_vector_layer(self, vector_layer, target_crs, output_path):
        """
        Reproject a vector layer to target CRS.

        Args:
            vector_layer: QgsVectorLayer to reproject
            target_crs: QgsCoordinateReferenceSystem target CRS
            output_path: Path for reprojected output

        Returns:
            QgsVectorLayer: Reprojected layer or None if failed
        """
        params_reproject = {
            'INPUT': vector_layer,
            'SOURCE_CRS': None,  # Auto-detect
            'TARGET_CRS': target_crs,
            'OPERATION': '',  # Use default
            'TARGET_EXTENT': None,
            'TARGET_EXTENT_CRS': None,
            'CONCAT_PREFIX': '',
            'NO_DATA': None,
            'COPY_SUBDATASETS': False,
            'DATA_TYPE': 0,
            'GPKG_LAYER_NAME': '',
            'OUTPUT': output_path
        }

        result = self.run_processing_algorithm(
            "native:reprojectlayer", params_reproject)
        if result and os.path.exists(output_path):
            reprojected_layer = QgsVectorLayer(
                output_path, f"{vector_layer.name()}_reprojected", "ogr")
            if reprojected_layer.isValid():
                self.log_message(
                    f"Vector layer reprojected successfully to {target_crs.authid()}")
                return reprojected_layer
            else:
                self.log_message(
                    "ERROR: Reprojected vector layer is not valid.")
                return None
        else:
            self.log_message("ERROR: Vector layer reprojection failed.")
            return None

    def process_layer_with_dem_clipping(self, layer, layer_name, output_filename, display=True):
        """
        Process a layer (vector or raster) by rasterizing if needed and clipping to DEM extent

        Args:
            layer: QgsMapLayer (vector or raster)
            layer_name: Name for logging
            output_filename: Output filename in geometry folder
            display: bool, whether to load the processed layer to QGIS (default: True)

        Returns:
            str: Path to processed layer or None if failed
        """
        if not layer:
            QMessageBox.warning(self.dialog, "Input Error",
                                f"Please select a {layer_name} layer.")
            return None

        if not self.check_prerequisites():
            return None

        # Get input CRS and DEM info
        input_crs = self.dialog.get_crs()

        # Use reprojected DEM layer if available, otherwise get from combo box
        if not self.dem_layer:
            dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
            # Check and reproject if needed
            self.check_and_reproject_dem(dem_layer)

        # Use the reprojected DEM layer stored in self.dem_layer
        dem_layer = self.dem_layer if self.dem_layer else self.dialog.mMapLayerComboBox_dem.currentLayer()

        # Get DEM extent and resolution from reprojected DEM
        extent_str, pixel_size_x, pixel_size_y = self.get_dem_extent_and_resolution()
        if not extent_str:
            return None

        # Get reprojected DEM extent with CRS for clipping
        dem_extent = self.dem_layer.extent()
        dem_crs = self.dem_layer.crs()
        # Format: xmin,xmax,ymin,ymax [EPSG:xxxxx]
        extent_projwin = f"{dem_extent.xMinimum()},{dem_extent.xMaximum()},{dem_extent.yMinimum()},{dem_extent.yMaximum()} [{dem_crs.authid()}]"

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        output_path = os.path.join(geometry_folder, output_filename)

        # Check if already processed
        if os.path.exists(output_path):
            self.log_message(f"{layer_name} already processed.")
            if display:
                self.load_layer(output_path, layer_name)
            return output_path

        self.log_message(f"Processing {layer_name}...")

        # Determine if layer is vector or raster
        if isinstance(layer, QgsVectorLayer):
            self.log_message(
                f"Vector layer detected. Checking CRS for {layer_name}...")

            # Check if vector layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(
                        f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(
                        geometry_folder, f"temp_{layer_name}_reprojected.gpkg")
                    reprojected_layer = self.reproject_vector_layer(
                        layer, input_crs, temp_reprojected_path)
                    if reprojected_layer:
                        layer = reprojected_layer
                    else:
                        self.log_message(
                            f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(
                        f"{layer_name} CRS matches input CRS. No reprojection needed.")

            # Rasterize vector layer using DEM's georeferenced units
            self.log_message(f"Rasterizing {layer_name}...")

            # Get the first field for rasterization (or use a default)
            fields = layer.fields()
            field_name = fields[0].name() if fields else 'id'

            # Calculate width and height from reprojected DEM pixel size
            dem_extent = self.dem_layer.extent()
            width = int((dem_extent.xMaximum() -
                        dem_extent.xMinimum()) / pixel_size_x)
            height = int((dem_extent.yMaximum() -
                         dem_extent.yMinimum()) / pixel_size_y)

            params_rasterize = {
                'INPUT': layer,
                'FIELD': field_name,
                'BURN': 0,
                'USE_Z': False,
                'UNITS': 1,  # Georeferenced units
                'WIDTH': width,
                'HEIGHT': height,
                'EXTENT': extent_str,
                'NODATA': -9999,
                'OPTIONS': None,
                'DATA_TYPE': 5,  # Float32
                'INIT': None,
                'INVERT': False,
                'EXTRA': '',
                'OUTPUT': output_path
            }

            result = self.run_processing_algorithm(
                "gdal:rasterize", params_rasterize)

        elif isinstance(layer, QgsRasterLayer):
            self.log_message(
                f"Raster layer detected. Clipping {layer_name} to DEM extent...")

            # Check if raster layer CRS matches input CRS
            layer_crs = layer.crs()
            if input_crs.isValid() and layer_crs.isValid():
                if layer_crs.authid() != input_crs.authid():
                    self.log_message(
                        f"{layer_name} CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                    temp_reprojected_path = os.path.join(
                        geometry_folder, f"temp_{layer_name}_reprojected.tif")

                    params_warp = {
                        'INPUT': layer.source(),
                        'SOURCE_CRS': None,
                        'TARGET_CRS': input_crs,
                        'RESAMPLING': 0,  # Nearest neighbor
                        'NODATA': -9999,
                        'TARGET_RESOLUTION': None,
                        'OPTIONS': None,
                        'DATA_TYPE': 0,
                        'TARGET_EXTENT': None,
                        'TARGET_EXTENT_CRS': None,
                        'MULTITHREADING': False,
                        'EXTRA': '',
                        'OUTPUT': temp_reprojected_path
                    }

                    warp_result = self.run_processing_algorithm(
                        "gdal:warpreproject", params_warp)
                    if warp_result and os.path.exists(temp_reprojected_path):
                        layer = QgsRasterLayer(
                            temp_reprojected_path, f"{layer_name}_reprojected")
                        self.log_message(
                            f"{layer_name} reprojected successfully.")
                    else:
                        self.log_message(
                            f"WARNING: Reprojection failed. Using original {layer_name} layer.")
                else:
                    self.log_message(
                        f"{layer_name} CRS matches input CRS. No reprojection needed.")

            # Clip raster to DEM extent using reprojected DEM extent
            params_clip = {
                'INPUT': layer,
                'PROJWIN': extent_projwin,  # Extent from reprojected DEM with CRS
                'OVERCRS': False,
                'NODATA': None,
                'OPTIONS': None,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'OUTPUT': output_path
            }

            result = self.run_processing_algorithm(
                "gdal:cliprasterbyextent", params_clip)

        else:
            self.log_message(f"ERROR: Unsupported layer type for {layer_name}")
            return None

        if result:
            if display:
                self.load_layer(output_path, layer_name)
            self.log_message(
                f"{layer_name} processing completed successfully.")
            return output_path
        else:
            self.log_message(f"{layer_name} processing failed.")
            return None

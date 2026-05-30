# -*- coding: utf-8 -*-
"""Geology class raster preparation."""
from ..common import (
    os,
    QMessageBox,
    QgsVectorLayer,
    QgsRasterLayer,
)


class GeologyProcessingMixin:
    """Geology class raster preparation."""

    def process_geology(self):
        """Process Geology layer by rasterizing using GDAL with filled DEM extent and resolution"""
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Check if filled DEM exists
        if not self._ensure_filled_dem():
            return
        
        # Get geology layer
        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology layer.")
            return
        
        # Check if it's a vector layer (rasterization is for vector layers)
        if not isinstance(layer, QgsVectorLayer):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology layer must be a vector layer for rasterization.")
            return
        
        # Check if output already exists
        geometry_folder = os.path.join(self.dialog.project_folder, "Geometry")
        self.geology_path = os.path.join(geometry_folder, "3_geology_processed.tif")
        
        if os.path.exists(self.geology_path):
            self.log_message("Geology layer already processed. Loading existing file...")
            self.load_layer(self.geology_path, "Geology")
            return
        
        self.log_message("Processing Geology layer...")
        
        # Get input CRS from dialog
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        # Check if geology layer CRS matches input CRS
        layer_crs = layer.crs()
        temp_reprojected_path = None
        layer_to_rasterize = layer
        
        if layer_crs.isValid():
            if layer_crs.authid() != input_crs.authid():
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                temp_reprojected_path = os.path.join(
                    geometry_folder, "temp_geology_reprojected.shp")
                
                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_rasterize = reprojected_layer
                    self.log_message(
                        f"Geology layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "WARNING: Reprojection failed. Using original geology layer.")
                    temp_reprojected_path = None
            else:
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Geology layer CRS is not valid. Proceeding with rasterization...")
        
        # Get filled DEM layer for extent and resolution
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    os.remove(temp_reprojected_path)
                except:
                    pass
            return
        
        # Get extent and dimensions from filled DEM
        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()  # Pixel width
        height = filled_dem_layer.height()  # Pixel height

        cell_size_width = (raster_extent.xMaximum() -
                           raster_extent.xMinimum()) / width
        cell_size_height = (raster_extent.yMaximum() -
                            raster_extent.yMinimum()) / height
        
        # Format extent as string: "xmin,xmax,ymin,ymax [EPSG:xxxxx]"
        extent_str = f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},{raster_extent.yMinimum()},{raster_extent.yMaximum()} [{input_crs.authid()}]"
        
        self.log_message(f"Using filled DEM extent: {extent_str}")
        self.log_message(f"Using filled DEM resolution: {width} x {height} pixels")
        
        # Get field name for rasterization (try common field names, or use first field)
        fields = layer_to_rasterize.fields()
        field_name = None
        # Try common geology field names
        for field_candidate in ['GEO_CLASS', 'geology', 'Geology', 'GEO', 'CLASS', 'class', 'id']:
            if field_candidate in [f.name() for f in fields]:
                field_name = field_candidate
                break
        
        # If no common field found, use first field
        if not field_name and fields:
            field_name = fields[0].name()
        
        if not field_name:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology layer must have at least one attribute field for rasterization.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    os.remove(temp_reprojected_path)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{field_name}' for rasterization")
        
        # Get layer source path (for vector layers)
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            # Use the reprojected layer path
            layer_input = temp_reprojected_path
        else:
            layer_source = layer_to_rasterize.source()
            if not layer_source or not os.path.exists(layer_source):
                # If source is not a file path, use the layer directly
                layer_input = layer_to_rasterize
            else:
                layer_input = layer_source
        
        # Rasterize using GDAL
        params_rasterize = {
            'INPUT': layer_input,
            'FIELD': field_name,
            'BURN': 0,
            'USE_Z': False,
            'UNITS': 1,  # Georeferenced units
            'WIDTH': cell_size_width,
            'HEIGHT': cell_size_height,
            'EXTENT': extent_str,   
            'NODATA': -9999,
            'OPTIONS': None,
            'DATA_TYPE': 5,  # Int32
            'INIT': None,
            'INVERT': False,
            'EXTRA': '',
            'OUTPUT': self.geology_path
        }
        
        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        if result:
            self.load_layer(self.geology_path, "Geology")
            self.log_message("Geology layer rasterized successfully.")
        else:
            self.log_message("ERROR: Geology rasterization failed.")
            self.geology_path = None
        
        # Clean up temporary reprojected file if it was created
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                # Remove shapefile and associated files (.shp, .shx, .dbf, .prj, etc.)
                base_name = os.path.splitext(temp_reprojected_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary reprojected geology file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary file: {e}")

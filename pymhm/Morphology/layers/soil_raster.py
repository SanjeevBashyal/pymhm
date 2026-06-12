# -*- coding: utf-8 -*-
"""Soil raster preparation and class rasterization."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsVectorLayer,
    QgsRasterLayer,
    QgsFeature,
    QgsFields,
    create_vector_file_writer,
    qgs_field,
    processing,
)
from ..watershed.dem_fill import DemFillMixin
from ..core.predecessors import PredecessorMixin
from ..core.layer_preparation import LayerPreparationMixin
from .soil_lookup import SoilLookupMixin


class SoilRasterMixin(
        LayerPreparationMixin,
        DemFillMixin,
        PredecessorMixin,
        SoilLookupMixin):
    """Soil raster preparation and class rasterization."""

    def process_soil(self) -> None:
        """Process Soil layer by rasterizing using GDAL with filled DEM extent and resolution"""
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Check if filled DEM exists
        if not self._ensure_filled_dem(self.fill_dem):
            return
        
        # Get soil layer
        layer = self.dialog.mMapLayerComboBox_soil.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil layer.")
            return
        
        # Check if it's a vector layer (rasterization is for vector layers)
        if not isinstance(layer, QgsVectorLayer):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Soil layer must be a vector layer for rasterization.")
            return
        
        # Check if output already exists
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        output_path = os.path.join(geometry_folder, "3_soil.tif")
        
        if os.path.exists(output_path):
            self.log_message("Soil layer already processed. Loading existing file...")
            self.load_layer(output_path, "3_Soil")
            return
        
        self.log_message("Processing Soil layer...")
        
        # Get input CRS from dialog
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return
        
        # Check if soil layer CRS matches input CRS
        layer_crs = layer.crs()
        temp_reprojected_path = None
        layer_to_process = layer
        
        if layer_crs.isValid():
            if layer_crs.authid() != input_crs.authid():
                self.log_message(
                    f"Soil layer CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                temp_reprojected_path = os.path.join(
                    geometry_folder, "temp_soil_reprojected.shp")
                
                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_process = reprojected_layer
                    self.log_message(
                        f"Soil layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "WARNING: Reprojection failed. Using original soil layer.")
                    temp_reprojected_path = None
            else:
                self.log_message(
                    f"Soil layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Soil layer CRS is not valid. Proceeding with processing...")
        
        # Load soil lookup table
        lookup_mapping = self.load_soil_lookup_table()
        if not lookup_mapping:
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # The selected lookup field must also exist on the soil GIS layer.
        fields = layer_to_process.fields()
        field_names = [f.name() for f in fields]
        lookup_key_field = getattr(self, "soil_lookup_key_field", None)
        soil_key_field = (
            self._required_lookup_field(field_names, lookup_key_field)
            if lookup_key_field else None
        )
        
        if not soil_key_field:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Soil layer must contain the selected lookup field '{lookup_key_field}'.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{soil_key_field}' for soil class lookup")
        
        # Create a temporary layer with CLASS field
        temp_soil_with_class_path = os.path.join(
            geometry_folder, "temp_soil_with_class.shp")
        
        # Add CLASS field and populate it
        self.log_message("Adding CLASS field to soil layer...")
        
        # Create output fields (copy all existing fields + add CLASS)
        output_fields = QgsFields()
        for field in fields:
            output_fields.append(field)
        output_fields.append(qgs_field("CLASS", "Int"))
        
        # Remove existing file if it exists
        if os.path.exists(temp_soil_with_class_path):
            try:
                base_name = os.path.splitext(temp_soil_with_class_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
            except:
                pass
        
        # Get geometry type from input layer
        geometry_type = layer_to_process.wkbType()
        
        writer = create_vector_file_writer(
            temp_soil_with_class_path,
            output_fields,
            geometry_type,
            input_crs,
        )
        
        if writer.hasError():
            self.log_message(
                f"ERROR creating temporary soil layer: {writer.errorMessage()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error creating temporary soil layer:\n{writer.errorMessage()}")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Process features and add CLASS value
        features_processed = 0
        features_mapped = 0
        features_unmapped = 0
        
        for feature in layer_to_process.getFeatures():
            new_feat = QgsFeature(output_fields)
            new_feat.setGeometry(feature.geometry())
            
            # Copy all attributes
            attrs = feature.attributes()
            for i, attr in enumerate(attrs):
                if i < len(output_fields) - 1:  # Exclude the last field (CLASS)
                    new_feat.setAttribute(i, attr)
            
            # Get selected lookup value and map it to CLASS.
            soil_key_value = feature.attribute(soil_key_field)
            class_code = None
            
            soil_key = self._normalise_lookup_key(soil_key_value)
            if soil_key:
                class_code = lookup_mapping.get(soil_key)
            
            if class_code is not None:
                new_feat.setAttribute("CLASS", class_code)
                features_mapped += 1
            else:
                self.log_message(
                    f"ERROR: No CLASS mapping found for {soil_key_field}='{soil_key_value}'.")
                new_feat.setAttribute("CLASS", 0)
                features_unmapped += 1
            
            writer.addFeature(new_feat)
            features_processed += 1
        
        del writer
        
        self.log_message(
            f"Processed {features_processed} features. "
            f"Mapped: {features_mapped}, Unmapped: {features_unmapped}")
        
        if features_mapped == 0 or features_unmapped > 0:
            QMessageBox.warning(
                self.dialog, "Mapping Error",
                "Soil mapping failed. Every soil feature must map to a CLASS value in the lookup table.")
            # Clean up temporary files
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    base_name = os.path.splitext(temp_soil_with_class_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            return
        
        # Get filled DEM layer for extent and resolution
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            # Clean up temporary files
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    base_name = os.path.splitext(temp_reprojected_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    base_name = os.path.splitext(temp_soil_with_class_path)[0]
                    for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                        file_to_remove = base_name + ext
                        if os.path.exists(file_to_remove):
                            os.remove(file_to_remove)
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
        
        # Rasterize using GDAL with CLASS field
        params_rasterize = {
            'INPUT': temp_soil_with_class_path,
            'FIELD': 'CLASS',
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
            'OUTPUT': output_path
        }
        
        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        if result:
            self.load_layer(output_path, "3_Soil")
            self.log_message("Soil layer rasterized successfully.")
        else:
            self.log_message("ERROR: Soil rasterization failed.")
        
        # Clean up temporary files
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                base_name = os.path.splitext(temp_reprojected_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary reprojected soil file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary reprojected file: {e}")
        
        if os.path.exists(temp_soil_with_class_path):
            try:
                base_name = os.path.splitext(temp_soil_with_class_path)[0]
                for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj']:
                    file_to_remove = base_name + ext
                    if os.path.exists(file_to_remove):
                        os.remove(file_to_remove)
                self.log_message("Cleaned up temporary soil with class file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary soil with class file: {e}")

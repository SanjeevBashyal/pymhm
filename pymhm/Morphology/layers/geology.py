# -*- coding: utf-8 -*-
"""Geology class raster preparation."""
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
)
from ..watershed.dem_fill import DemFillMixin
from ..core.predecessors import PredecessorMixin
from ..core.layer_preparation import LayerPreparationMixin
from .lookup_fields import LookupFieldMixin


class GeologyProcessingMixin(
        LayerPreparationMixin,
        DemFillMixin,
        PredecessorMixin,
        LookupFieldMixin):
    """Geology class raster preparation."""

    def process_geology(self) -> None:
        """Process Geology layer by rasterizing using GDAL with filled DEM extent and resolution"""
        # Check prerequisites
        if not self.check_prerequisites():
            return
        
        # Check if filled DEM exists
        if not self._ensure_filled_dem(self.fill_dem):
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
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
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
        temp_geology_with_class_path = None
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
                    self._remove_geology_vector_dataset(temp_reprojected_path)
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
        
        # The selected lookup field must exist on both the lookup table and geology layer.
        fields = layer_to_rasterize.fields()
        field_names = [f.name() for f in fields]
        geology_lookup_mapping, geology_lookup_key_field = (
            self._load_geology_lookup_mapping())
        if not geology_lookup_mapping or not geology_lookup_key_field:
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                self._remove_geology_vector_dataset(temp_reprojected_path)
            return

        geology_key_field = self._required_lookup_field(
            field_names, geology_lookup_key_field)
        
        if not geology_key_field:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Geology layer must contain the selected lookup field '{geology_lookup_key_field}'.")
            # Clean up temporary file if it was created
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    self._remove_geology_vector_dataset(temp_reprojected_path)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{geology_key_field}' for geology class lookup")

        temp_geology_with_class_path = os.path.join(
            geometry_folder, "temp_geology_with_class.shp")
        layer_input = self._write_geology_class_layer(
            layer_to_rasterize,
            fields,
            geology_key_field,
            geology_lookup_mapping,
            temp_geology_with_class_path,
            input_crs,
        )
        if not layer_input:
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                self._remove_geology_vector_dataset(temp_reprojected_path)
            return
        
        # Rasterize using GDAL
        params_rasterize = {
            'INPUT': layer_input,
            'FIELD': 'MHM_CLASS',
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
                self._remove_geology_vector_dataset(temp_reprojected_path)
                self.log_message("Cleaned up temporary reprojected geology file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary file: {e}")

        if temp_geology_with_class_path and os.path.exists(temp_geology_with_class_path):
            try:
                self._remove_geology_vector_dataset(temp_geology_with_class_path)
                self.log_message("Cleaned up temporary geology with class file.")
            except Exception as e:
                self.log_message(
                    f"WARNING: Could not delete temporary geology with class file: {e}")

    def _load_geology_lookup_mapping(self):
        """Return geology lookup mapping and selected table key field."""
        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_geologyLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )
        if not lookup_layer:
            return None, None

        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Geology lookup layer is not valid.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a valid Geology lookup table layer.")
            return None, None

        field_names = lookup_layer.fields().names()
        key_field = self._selected_lookup_field(
            "comboBox_geologyLookupField", field_names)
        if not key_field:
            self.log_message(
                "ERROR: Please select a geology lookup field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select the geology lookup field that maps the geology layer to the lookup table.")
            return None, None

        class_field = self._required_lookup_field(field_names, "Geo-Class")
        if not class_field:
            self.log_message("ERROR: Geology lookup layer must contain a Geo-Class field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology lookup layer must contain a Geo-Class field.")
            return None, None

        self.log_message(
            f"Using geology lookup fields '{key_field}' -> '{class_field}'")
        lookup_mapping = {}
        for feature_index, feature in enumerate(lookup_layer.getFeatures(), start=1):
            try:
                lookup_key = self._normalise_lookup_key(feature[key_field])
                if not lookup_key:
                    self.log_message(
                        f"WARNING: Geology lookup row {feature_index} has an empty key. Skipping.")
                    continue

                class_code = self._coerce_lookup_int(feature[class_field])
                if class_code is None or class_code <= 0:
                    self.log_message(
                        f"WARNING: Geology lookup row {feature_index} has invalid Geo-Class '{feature[class_field]}'. Skipping.")
                    continue

                lookup_mapping[lookup_key] = class_code
            except Exception as e:
                self.log_message(
                    f"WARNING: Error parsing geology lookup feature {feature_index}: {e}")

        if lookup_mapping:
            self.log_message(
                f"Loaded {len(lookup_mapping)} entries from geology lookup layer.")
            return lookup_mapping, key_field

        self.log_message("ERROR: Geology lookup layer is empty or could not be parsed.")
        QMessageBox.warning(
            self.dialog, "Input Error",
            "Geology lookup layer is empty or could not be parsed.")
        return None, None

    def _write_geology_class_layer(
            self,
            source_layer,
            source_fields,
            source_field_name,
            lookup_mapping,
            output_path,
            output_crs):
        """Create a temporary geology vector with a numeric MHM_CLASS field."""
        self.log_message("Adding MHM_CLASS field to geology layer from lookup table...")
        self._remove_geology_vector_dataset(output_path)

        output_fields = QgsFields()
        for field in source_fields:
            output_fields.append(field)
        output_fields.append(qgs_field("MHM_CLASS", "Int"))

        writer = create_vector_file_writer(
            output_path,
            output_fields,
            source_layer.wkbType(),
            output_crs,
        )
        if writer.hasError():
            self.log_message(
                f"ERROR creating temporary geology layer: {writer.errorMessage()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error creating temporary geology layer:\n{writer.errorMessage()}")
            return None

        features_processed = 0
        features_mapped = 0
        features_unmapped = 0

        for feature in source_layer.getFeatures():
            new_feature = QgsFeature(output_fields)
            new_feature.setGeometry(feature.geometry())

            attrs = feature.attributes()
            for i, attr in enumerate(attrs):
                if i < len(output_fields) - 1:
                    new_feature.setAttribute(i, attr)

            lookup_key = self._normalise_lookup_key(
                feature.attribute(source_field_name))
            class_code = lookup_mapping.get(lookup_key) if lookup_key else None
            if class_code is not None:
                new_feature.setAttribute("MHM_CLASS", class_code)
                features_mapped += 1
            else:
                self.log_message(
                    f"ERROR: No Geo-Class mapping found for {source_field_name}='{feature.attribute(source_field_name)}'.")
                new_feature.setAttribute("MHM_CLASS", 0)
                features_unmapped += 1

            writer.addFeature(new_feature)
            features_processed += 1

        del writer

        self.log_message(
            f"Processed {features_processed} geology features. "
            f"Mapped: {features_mapped}, Unmapped: {features_unmapped}")
        if features_mapped == 0 or features_unmapped > 0:
            QMessageBox.warning(
                self.dialog, "Mapping Error",
                "Geology mapping failed. Every geology feature must map to a Geo-Class value in the lookup table.")
            self._remove_geology_vector_dataset(output_path)
            return None

        return output_path

    def _remove_geology_vector_dataset(self, path):
        """Remove a temporary shapefile dataset and common sidecar files."""
        base_name = os.path.splitext(path)[0]
        for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj', '.fix']:
            file_to_remove = base_name + ext
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)

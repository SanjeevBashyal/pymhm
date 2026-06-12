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

    def process_geology(self, write_classdefinition=True) -> None:
        """Process geology input into an mHM geology class raster."""
        if not self.check_prerequisites():
            return

        if not self._ensure_filled_dem(self.fill_dem):
            return

        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology layer.")
            return

        if not isinstance(layer, (QgsVectorLayer, QgsRasterLayer)):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology layer must be a vector or raster layer.")
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        os.makedirs(geometry_folder, exist_ok=True)
        self.geology_path = os.path.join(
            geometry_folder, "3_geology_processed.tif")

        if os.path.exists(self.geology_path):
            self.log_message("Geology layer already processed. Loading existing file...")
            self.load_layer(self.geology_path, "Geology")
            if write_classdefinition:
                self._write_geology_classdefinition_after_processing()
            return

        self.log_message("Processing Geology layer...")

        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return

        geology_lookup_mapping, geology_lookup_key_field = (
            self._load_geology_lookup_mapping())
        if not geology_lookup_mapping or not geology_lookup_key_field:
            return

        if isinstance(layer, QgsRasterLayer):
            geology_raster_created = self._process_geology_raster_layer(
                layer,
                geology_lookup_mapping,
                self.geology_path,
                input_crs,
                geometry_folder,
            )
            if geology_raster_created and write_classdefinition:
                self._write_geology_classdefinition_after_processing()
            return

        layer_crs = layer.crs()
        temp_reprojected_path = None
        layer_to_rasterize = layer

        if layer_crs.isValid():
            if layer_crs.authid() != input_crs.authid():
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) differs from input CRS ({input_crs.authid()}). Reprojecting...")
                temp_reprojected_path = os.path.join(
                    geometry_folder, "temp_geology_reprojected.gpkg")

                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_rasterize = reprojected_layer
                    self.log_message(
                        f"Geology layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "ERROR: Geology reprojection failed. Aborting geology raster preparation.")
                    QMessageBox.critical(
                        self.dialog,
                        "Reprojection Error",
                        "Geology layer reprojection failed. The geology raster was not prepared.")
                    return
            else:
                self.log_message(
                    f"Geology layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Geology layer CRS is not valid. Proceeding with rasterization...")

        fields = layer_to_rasterize.fields()
        field_names = [field.name() for field in fields]
        geology_key_field = self._required_lookup_field(
            field_names, geology_lookup_key_field)

        if not geology_key_field:
            QMessageBox.warning(
                self.dialog, "Input Error",
                f"Geology layer must contain the selected lookup field '{geology_lookup_key_field}'.")
            self._cleanup_geology_vector(temp_reprojected_path)
            return

        self.log_message(
            f"Using field '{geology_key_field}' for geology class lookup")

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
            self._cleanup_geology_vector(temp_reprojected_path)
            return

        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            self._cleanup_geology_vector(temp_reprojected_path)
            self._cleanup_geology_vector(temp_geology_with_class_path)
            return

        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        cell_size_width = (
            raster_extent.xMaximum() - raster_extent.xMinimum()) / width
        cell_size_height = (
            raster_extent.yMaximum() - raster_extent.yMinimum()) / height
        extent_str = (
            f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},"
            f"{raster_extent.yMinimum()},{raster_extent.yMaximum()} "
            f"[{input_crs.authid()}]"
        )

        self.log_message(f"Using filled DEM extent: {extent_str}")
        self.log_message(f"Using filled DEM resolution: {width} x {height} pixels")

        params_rasterize = {
            'INPUT': layer_input,
            'FIELD': 'MHM_CLASS',
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
            'OUTPUT': self.geology_path
        }

        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        geology_raster_created = result and os.path.exists(self.geology_path)
        if geology_raster_created:
            self.load_layer(self.geology_path, "Geology")
            self.log_message("Geology layer rasterized successfully.")
        else:
            self.log_message("ERROR: Geology rasterization failed.")
            self.geology_path = None

        self._cleanup_geology_vector(temp_reprojected_path)
        self._cleanup_geology_vector(temp_geology_with_class_path)

        if geology_raster_created and write_classdefinition:
            self._write_geology_classdefinition_after_processing()

    def _load_geology_lookup_mapping(self):
        """Return geology lookup mapping and selected table key field."""
        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_geologyLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )
        if not lookup_layer:
            self.log_message("ERROR: Geology lookup table layer not selected.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Geology lookup table layer.")
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
            self.log_message("ERROR: Please select a geology lookup field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select the geology lookup field that maps the geology layer to the lookup table.")
            return None, None

        class_field = self._required_lookup_field(field_names, "GEOLOGY_CLASS")
        if not class_field:
            self.log_message(
                "ERROR: Geology lookup layer must contain a GEOLOGY_CLASS field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Geology lookup layer must contain a GEOLOGY_CLASS field.")
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
                        f"WARNING: Geology lookup row {feature_index} has invalid GEOLOGY_CLASS '{feature[class_field]}'. Skipping.")
                    continue

                lookup_mapping[lookup_key] = class_code
                if len(lookup_mapping) <= 10:
                    self.log_message(f"  Mapped: {lookup_key} -> {class_code}")
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
        self.log_message("Adding mapped geology class field to geology layer...")
        self._remove_vector_output(output_path)

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
            for index, attr in enumerate(attrs):
                if index < len(output_fields) - 1:
                    new_feature.setAttribute(index, attr)

            lookup_key = self._normalise_lookup_key(
                feature.attribute(source_field_name))
            class_code = lookup_mapping.get(lookup_key) if lookup_key else None
            if class_code is not None:
                new_feature.setAttribute("MHM_CLASS", class_code)
                features_mapped += 1
            else:
                self.log_message(
                    f"ERROR: No GEOLOGY_CLASS mapping found for {source_field_name}='{feature.attribute(source_field_name)}'.")
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
                "Geology mapping failed. Every geology feature must map to a GEOLOGY_CLASS value in the lookup table.")
            self._remove_vector_output(output_path)
            return None

        return output_path

    def _process_geology_raster_layer(
            self,
            layer,
            lookup_mapping,
            output_path,
            input_crs,
            geometry_folder):
        """Prepare a geology class raster from a categorical raster input."""
        self.log_message("Processing raster geology layer with lookup table...")

        reclass_table = []
        for lookup_key, class_code in sorted(
                lookup_mapping.items(),
                key=lambda item: self._numeric_sort_key(item[0])):
            try:
                raster_value = float(lookup_key)
            except (TypeError, ValueError):
                self.log_message(
                    "ERROR: Raster geology inputs require numeric values in the selected lookup field. "
                    f"Could not convert '{lookup_key}'.")
                QMessageBox.warning(
                    self.dialog,
                    "Input Error",
                    "Raster geology inputs require numeric values in the selected lookup field.")
                return False

            reclass_table.extend([
                str(raster_value - 0.5),
                str(raster_value + 0.5),
                str(int(class_code)),
            ])

        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            self.log_message("ERROR: Cannot read filled DEM layer.")
            QMessageBox.critical(
                self.dialog, "Error",
                "Cannot read filled DEM layer.")
            return False

        raster_extent = filled_dem_layer.extent()
        width = filled_dem_layer.width()
        height = filled_dem_layer.height()
        cell_size_width = (
            raster_extent.xMaximum() - raster_extent.xMinimum()) / width
        cell_size_height = (
            raster_extent.yMaximum() - raster_extent.yMinimum()) / height
        extent_str = (
            f"{raster_extent.xMinimum()},{raster_extent.xMaximum()},"
            f"{raster_extent.yMinimum()},{raster_extent.yMaximum()} "
            f"[{input_crs.authid()}]"
        )

        aligned_path = os.path.join(
            geometry_folder, "temp_geology_raster_aligned.tif")
        params_warp = {
            'INPUT': layer.source(),
            'SOURCE_CRS': None,
            'TARGET_CRS': input_crs,
            'RESAMPLING': 0,
            'NODATA': -9999,
            'TARGET_RESOLUTION': min(cell_size_width, cell_size_height),
            'OPTIONS': None,
            'DATA_TYPE': 0,
            'TARGET_EXTENT': extent_str,
            'TARGET_EXTENT_CRS': input_crs,
            'MULTITHREADING': False,
            'EXTRA': '',
            'OUTPUT': aligned_path
        }

        warp_result = self.run_processing_algorithm(
            "gdal:warpreproject", params_warp)
        if not warp_result or not os.path.exists(aligned_path):
            self.log_message("ERROR: Failed to align geology raster to filled DEM grid.")
            return False

        params_reclass = {
            'INPUT_RASTER': aligned_path,
            'RASTER_BAND': 1,
            'TABLE': reclass_table,
            'NO_DATA': -9999,
            'RANGE_BOUNDARIES': 0,
            'NODATA_FOR_MISSING': True,
            'DATA_TYPE': 5,
            'CREATE_OPTIONS': None,
            'OUTPUT': output_path
        }

        result = self.run_processing_algorithm(
            "native:reclassifybytable", params_reclass)
        geology_raster_created = result and os.path.exists(output_path)
        if geology_raster_created:
            self.load_layer(output_path, "Geology")
            self.log_message("Raster geology layer reclassified successfully.")
        else:
            self.log_message("ERROR: Geology raster reclassification failed.")
            self.geology_path = None

        if os.path.exists(aligned_path):
            try:
                os.remove(aligned_path)
            except Exception as e:
                self.log_message(
                    f"WARNING: Could not delete temporary aligned geology raster: {e}")

        return bool(geology_raster_created)

    def _write_geology_classdefinition_after_processing(self):
        """Write the geology classdefinition when this processor includes the writer mixin."""
        writer = getattr(self, "geology_classification_writer", None)
        if writer is None:
            return None
        return writer()

    def _cleanup_geology_vector(self, path):
        """Remove a temporary geology vector dataset if it exists."""
        if path and os.path.exists(path):
            try:
                self._remove_vector_output(path)
            except Exception as e:
                self.log_message(
                    f"WARNING: Could not delete temporary geology vector '{path}': {e}")

    def _numeric_sort_key(self, value):
        """Sort lookup keys numerically when possible."""
        try:
            return (0, float(value))
        except (TypeError, ValueError):
            return (1, str(value))

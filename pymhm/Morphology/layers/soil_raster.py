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

    def process_soil(self, write_classdefinition=True) -> None:
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
        
        if not isinstance(layer, (QgsVectorLayer, QgsRasterLayer)):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Soil layer must be a vector or raster layer.")
            return
        
        # Check if output already exists
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        os.makedirs(geometry_folder, exist_ok=True)
        output_path = os.path.join(geometry_folder, "3_soil.tif")
        
        if os.path.exists(output_path):
            self.log_message("Soil layer already processed. Loading existing file...")
            self.load_layer(output_path, "3_Soil")
            if write_classdefinition:
                self._write_soil_classdefinition_after_processing()
            return
        
        self.log_message("Processing Soil layer...")
        
        # Get input CRS from dialog
        input_crs = self.dialog.get_crs()
        if not input_crs.isValid():
            QMessageBox.warning(
                self.dialog, "CRS Error",
                "Please set a valid input CRS.")
            return

        # Load soil lookup table
        lookup_mapping = self.load_soil_lookup_table()
        if not lookup_mapping:
            return

        if isinstance(layer, QgsRasterLayer):
            soil_raster_created = self._process_soil_raster_layer(
                layer,
                lookup_mapping,
                output_path,
                input_crs,
                geometry_folder,
            )
            if soil_raster_created and write_classdefinition:
                self._write_soil_classdefinition_after_processing()
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
                    geometry_folder, "temp_soil_reprojected.gpkg")
                
                reprojected_layer = self.reproject_vector_layer(
                    layer, input_crs, temp_reprojected_path)
                if reprojected_layer:
                    layer_to_process = reprojected_layer
                    self.log_message(
                        f"Soil layer reprojected successfully to {input_crs.authid()}")
                else:
                    self.log_message(
                        "ERROR: Soil reprojection failed. Aborting soil raster preparation.")
                    QMessageBox.critical(
                        self.dialog,
                        "Reprojection Error",
                        "Soil layer reprojection failed. The soil raster was not prepared.")
                    temp_reprojected_path = None
                    return
            else:
                self.log_message(
                    f"Soil layer CRS ({layer_crs.authid()}) matches input CRS. No reprojection needed.")
        else:
            self.log_message(
                "WARNING: Soil layer CRS is not valid. Proceeding with processing...")
        
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
                    self._remove_vector_output(temp_reprojected_path)
                except:
                    pass
            return
        
        self.log_message(f"Using field '{soil_key_field}' for soil class lookup")
        
        # Create a temporary layer with mapped mHM soil class field.
        temp_soil_with_class_path = os.path.join(
            geometry_folder, "temp_soil_with_class.shp")
        
        # Add internal class field and populate it
        self.log_message("Adding mapped soil class field to soil layer...")
        
        # Create output fields (copy all existing fields + add internal class)
        output_fields = QgsFields()
        for field in fields:
            output_fields.append(field)
        output_fields.append(qgs_field("MHM_CLASS", "Int"))
        
        self._remove_vector_output(temp_soil_with_class_path)
        
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
                    self._remove_vector_output(temp_reprojected_path)
                except:
                    pass
            return
        
        # Process features and add mapped soil class value
        features_processed = 0
        features_mapped = 0
        features_unmapped = 0
        
        for feature in layer_to_process.getFeatures():
            new_feat = QgsFeature(output_fields)
            new_feat.setGeometry(feature.geometry())
            
            # Copy all attributes
            attrs = feature.attributes()
            for i, attr in enumerate(attrs):
                if i < len(output_fields) - 1:  # Exclude internal class field
                    new_feat.setAttribute(i, attr)
            
            # Get selected lookup value and map it to SOIL_CLASS.
            soil_key_value = feature.attribute(soil_key_field)
            class_code = None
            
            soil_key = self._normalise_lookup_key(soil_key_value)
            if soil_key:
                class_code = lookup_mapping.get(soil_key)
            
            if class_code is not None:
                new_feat.setAttribute("MHM_CLASS", class_code)
                features_mapped += 1
            else:
                self.log_message(
                    f"ERROR: No SOIL_CLASS mapping found for {soil_key_field}='{soil_key_value}'.")
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
                "Soil mapping failed. Every soil feature must map to a SOIL_CLASS value in the lookup table.")
            # Clean up temporary files
            if temp_reprojected_path and os.path.exists(temp_reprojected_path):
                try:
                    self._remove_vector_output(temp_reprojected_path)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    self._remove_vector_output(temp_soil_with_class_path)
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
                    self._remove_vector_output(temp_reprojected_path)
                except:
                    pass
            if os.path.exists(temp_soil_with_class_path):
                try:
                    self._remove_vector_output(temp_soil_with_class_path)
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
        
        # Rasterize using GDAL with mapped class field
        params_rasterize = {
            'INPUT': temp_soil_with_class_path,
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
            'OUTPUT': output_path
        }
        
        result = self.run_processing_algorithm("gdal:rasterize", params_rasterize)
        soil_raster_created = result and os.path.exists(output_path)
        if soil_raster_created:
            self.load_layer(output_path, "3_Soil")
            self.log_message("Soil layer rasterized successfully.")
        else:
            self.log_message("ERROR: Soil rasterization failed.")
        
        # Clean up temporary files
        if temp_reprojected_path and os.path.exists(temp_reprojected_path):
            try:
                self._remove_vector_output(temp_reprojected_path)
                self.log_message("Cleaned up temporary reprojected soil file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary reprojected file: {e}")
        
        if os.path.exists(temp_soil_with_class_path):
            try:
                self._remove_vector_output(temp_soil_with_class_path)
                self.log_message("Cleaned up temporary soil with class file.")
            except Exception as e:
                self.log_message(f"WARNING: Could not delete temporary soil with class file: {e}")

        if soil_raster_created and write_classdefinition:
            self._write_soil_classdefinition_after_processing()

    def _process_soil_raster_layer(
            self,
            layer,
            lookup_mapping,
            output_path,
            input_crs,
            geometry_folder):
        """Prepare a soil class raster from a categorical raster input."""
        self.log_message("Processing raster soil layer with lookup table...")

        reclass_table = []
        for lookup_key, class_code in sorted(
                lookup_mapping.items(),
                key=lambda item: self._numeric_sort_key(item[0])):
            try:
                raster_value = float(lookup_key)
            except (TypeError, ValueError):
                self.log_message(
                    "ERROR: Raster soil inputs require numeric values in the selected lookup field. "
                    f"Could not convert '{lookup_key}'.")
                QMessageBox.warning(
                    self.dialog,
                    "Input Error",
                    "Raster soil inputs require numeric values in the selected lookup field.")
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

        aligned_path = os.path.join(geometry_folder, "temp_soil_raster_aligned.tif")
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
            self.log_message("ERROR: Failed to align soil raster to filled DEM grid.")
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
        soil_raster_created = result and os.path.exists(output_path)
        if soil_raster_created:
            self.load_layer(output_path, "3_Soil")
            self.log_message("Raster soil layer reclassified successfully.")
        else:
            self.log_message("ERROR: Soil raster reclassification failed.")

        if os.path.exists(aligned_path):
            try:
                os.remove(aligned_path)
            except Exception as e:
                self.log_message(
                    f"WARNING: Could not delete temporary aligned soil raster: {e}")

        return bool(soil_raster_created)

    def _write_soil_classdefinition_after_processing(self):
        """Write the soil classdefinition when this processor includes the writer mixin."""
        writer = getattr(self, "soil_classdefinition_writer", None)
        if writer is None:
            return None
        return writer()

    def _numeric_sort_key(self, value):
        """Sort lookup keys numerically when possible."""
        try:
            return (0, float(value))
        except (TypeError, ValueError):
            return (1, str(value))

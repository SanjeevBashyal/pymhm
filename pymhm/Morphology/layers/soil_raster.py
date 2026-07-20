# -*- coding: utf-8 -*-
"""Soil raster preparation through :mod:`mhm_tools`."""
from __future__ import annotations

from ..common import (
    QMessageBox,
    QgsRasterLayer,
    QgsVectorLayer,
    morph_folder,
    os,
    project_geometry_folder,
)
from ..core.layer_preparation import LayerPreparationMixin
from ..core.predecessors import PredecessorMixin
from ..core.soil_sources import (
    local_layer_source,
    materialize_vector_layer,
    remove_vector_dataset,
)
from ..watershed.dem_fill import DemFillMixin
from ...mhm_tools_to_integrate.setup_creation import (
    format_soil_file,
    rasterize_soil_map,
    write_soil_classdefinition_file,
)
from .soil_lookup import SoilLookupMixin


class SoilRasterMixin(
        LayerPreparationMixin,
        DemFillMixin,
        PredecessorMixin,
        SoilLookupMixin):
    """Prepare an mHM soil class raster from a QGIS vector or raster layer."""

    def process_soil(self, write_classdefinition=True) -> bool:
        """Create ``3_soil.tif`` on the filled-DEM grid."""
        if not self.check_prerequisites():
            return False
        if not self._ensure_filled_dem(self.fill_dem):
            return False

        layer = self.dialog.mMapLayerComboBox_soil.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error", "Please select a Soil layer.")
            return False
        if not isinstance(layer, (QgsVectorLayer, QgsRasterLayer)):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Soil layer must be a vector or raster layer.",
            )
            return False

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        static_morph_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(geometry_folder, exist_ok=True)
        os.makedirs(static_morph_folder, exist_ok=True)
        output_path = os.path.join(geometry_folder, "3_soil.tif")
        classdefinition_path = os.path.join(
            static_morph_folder, "soil_classdefinition.txt")

        if os.path.exists(output_path):
            self.log_message(
                "Soil layer already processed. Loading existing file...")
            self.load_layer(output_path, "3_Soil")
            if not write_classdefinition:
                return True
            return bool(self._write_soil_classdefinition_after_processing())

        lookup_mapping = self.load_soil_lookup_table()
        lookup_mapping_field = getattr(
            self, "soil_lookup_key_field", None)
        if not lookup_mapping or not lookup_mapping_field:
            return False

        lookup_combo = getattr(
            self.dialog, "mMapLayerComboBox_soilLookup", None)
        lookup_layer = (
            lookup_combo.currentLayer() if lookup_combo is not None else None)
        if not lookup_layer:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a Soil lookup table layer.",
            )
            return False

        soil_mapping_field = lookup_mapping_field
        if isinstance(layer, QgsVectorLayer):
            soil_mapping_field = self._required_lookup_field(
                layer.fields().names(), lookup_mapping_field)
            if not soil_mapping_field:
                QMessageBox.warning(
                    self.dialog,
                    "Input Error",
                    "Soil layer must contain the selected lookup field "
                    f"'{lookup_mapping_field}'.",
                )
                return False

        temp_lookup = os.path.join(geometry_folder, "temp_soil_lookup.gpkg")
        temp_vector = os.path.join(geometry_folder, "temp_soil_input.gpkg")
        temp_raster = os.path.join(geometry_folder, "temp_soil_input.tif")
        temp_classdefinition = os.path.join(
            geometry_folder, "temp_soil_classdefinition.txt")
        dem_crs = self._layer_crs_text(
            QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        )

        self.log_message("Processing Soil layer through mhm-tools...")
        try:
            lookup_path = local_layer_source(lookup_layer)
            if lookup_path is None:
                lookup_path = materialize_vector_layer(lookup_layer, temp_lookup)

            if isinstance(layer, QgsVectorLayer):
                input_path = local_layer_source(layer)
                if input_path is None:
                    input_path = materialize_vector_layer(layer, temp_vector)
                rasterize_soil_map(
                    input_file=input_path,
                    dem_file=self.filled_dem_path,
                    output_file=output_path,
                    mapping_field=soil_mapping_field,
                    lookup_table=lookup_path,
                    lookup_mapping_field=lookup_mapping_field,
                    input_crs=self._layer_crs_text(layer),
                    dem_crs=dem_crs,
                    log=self.log_message,
                )
                if write_classdefinition:
                    write_soil_classdefinition_file(
                        lookup_table=lookup_path,
                        output_file=temp_classdefinition,
                        log=self.log_message,
                    )
            else:
                input_path = local_layer_source(layer)
                if input_path is None:
                    input_path = self._materialize_soil_raster(
                        layer, temp_raster)
                format_soil_file(
                    input_file=input_path,
                    dem_file=self.filled_dem_path,
                    output_file=output_path,
                    lookup_table=lookup_path,
                    mapping_field=lookup_mapping_field,
                    classdefinition_file=temp_classdefinition,
                    input_crs=self._layer_crs_text(layer),
                    dem_crs=dem_crs,
                    log=self.log_message,
                )

            if not os.path.exists(output_path):
                raise RuntimeError("mhm-tools did not create the soil raster.")
            if write_classdefinition:
                if not os.path.exists(temp_classdefinition):
                    raise RuntimeError(
                        "mhm-tools did not create soil_classdefinition.txt.")
                os.replace(temp_classdefinition, classdefinition_path)
        except Exception as error:
            self.log_message(f"ERROR preparing soil data: {error}")
            QMessageBox.critical(
                self.dialog,
                "Soil Processing Error",
                f"Could not prepare the soil data:\n{error}",
            )
            if os.path.exists(output_path):
                os.remove(output_path)
            return False
        finally:
            self._cleanup_soil_temporary_inputs(
                temp_lookup,
                temp_vector,
                temp_raster,
                temp_classdefinition,
            )

        self.load_layer(output_path, "3_Soil")
        if write_classdefinition:
            self.mark_output_prepared(
                classdefinition_path,
                name="soil_classdefinition.txt",
                loaded=False,
            )
        self.log_message("Soil layer prepared successfully.")
        return True

    def _materialize_soil_raster(self, layer, output_path):
        """Write a QGIS raster, including provider URIs, to a GeoTIFF."""
        self._remove_soil_raster_output(output_path)
        result = self.run_processing_algorithm(
            "gdal:translate",
            {
                "INPUT": layer,
                "TARGET_CRS": None,
                "NODATA": None,
                "COPY_SUBDATASETS": False,
                "OPTIONS": None,
                "EXTRA": "",
                "DATA_TYPE": 0,
                "OUTPUT": output_path,
            },
        )
        if not result or not os.path.exists(output_path):
            raise RuntimeError("Could not materialize the selected soil raster.")
        return output_path

    def _cleanup_soil_temporary_inputs(
            self, lookup_path, vector_path, raster_path, definition_path):
        """Remove path-based copies made from the selected QGIS layers."""
        for vector_output in (lookup_path, vector_path):
            try:
                remove_vector_dataset(vector_output)
            except Exception as error:
                self.log_message(
                    "WARNING: Could not remove temporary soil vector "
                    f"'{vector_output}': {error}")
        try:
            self._remove_soil_raster_output(raster_path)
        except Exception as error:
            self.log_message(
                "WARNING: Could not remove temporary soil raster "
                f"'{raster_path}': {error}")
        try:
            if os.path.exists(definition_path):
                os.remove(definition_path)
        except Exception as error:
            self.log_message(
                "WARNING: Could not remove temporary soil classdefinition "
                f"'{definition_path}': {error}")

    def _remove_soil_raster_output(self, path):
        """Remove a raster and GDAL sidecars if they exist."""
        for candidate in (path, f"{path}.aux.xml", f"{path}.ovr"):
            if os.path.exists(candidate):
                os.remove(candidate)

    def _write_soil_classdefinition_after_processing(self):
        """Regenerate the definition for a cached soil raster."""
        writer = getattr(self, "soil_classdefinition_writer", None)
        if writer is None:
            self.log_message(
                "ERROR: Soil classdefinition writer is not available.")
            return None
        return writer()

    @staticmethod
    def _layer_crs_text(layer):
        """Return a CRS string suitable for assigning missing file metadata."""
        crs = layer.crs() if layer is not None else None
        if crs is None or not crs.isValid():
            return None
        authid = crs.authid()
        return authid or crs.toWkt()

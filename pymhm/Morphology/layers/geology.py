# -*- coding: utf-8 -*-
"""Geology raster preparation through :mod:`mhm_tools`."""
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
    format_geology_file,
    rasterize_geology_map,
    write_geology_classdefinition_file,
)
from .lookup_fields import LookupFieldMixin


class GeologyProcessingMixin(
    LayerPreparationMixin,
    DemFillMixin,
    PredecessorMixin,
    LookupFieldMixin,
):
    """Prepare an mHM geology-class raster from a vector or raster layer."""

    def process_geology(self, write_classdefinition=True) -> bool:
        """Create ``3_geology_processed.tif`` on the filled-DEM grid."""
        if not self.check_prerequisites():
            return False
        if not self._ensure_filled_dem(self.fill_dem):
            return False

        layer = self.dialog.mMapLayerComboBox_geology.currentLayer()
        if not layer:
            QMessageBox.warning(
                self.dialog, "Input Error", "Please select a Geology layer."
            )
            return False
        if not isinstance(layer, (QgsVectorLayer, QgsRasterLayer)):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Geology layer must be a vector or raster layer.",
            )
            return False

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        static_morph_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(geometry_folder, exist_ok=True)
        os.makedirs(static_morph_folder, exist_ok=True)
        output_path = os.path.join(geometry_folder, "3_geology_processed.tif")
        classdefinition_path = os.path.join(
            static_morph_folder, "geology_classdefinition.txt"
        )
        metadata_path = os.path.join(geometry_folder, "geology_class_metadata.json")
        self.geology_path = output_path

        if os.path.exists(output_path):
            self.log_message(
                "Geology layer already processed. Loading existing file..."
            )
            self.load_layer(output_path, "Geology")
            if write_classdefinition:
                return bool(self._write_geology_classdefinition_after_processing())
            return True

        lookup_layer, lookup_mapping_field = self._geology_lookup_input()
        if lookup_layer is None:
            return False

        geology_mapping_field = lookup_mapping_field
        if isinstance(layer, QgsVectorLayer):
            geology_mapping_field = self._required_lookup_field(
                layer.fields().names(), lookup_mapping_field
            )
            if not geology_mapping_field:
                QMessageBox.warning(
                    self.dialog,
                    "Input Error",
                    "Geology layer must contain the selected lookup field "
                    f"'{lookup_mapping_field}'.",
                )
                return False

        temp_lookup = os.path.join(geometry_folder, "temp_geology_lookup.gpkg")
        temp_vector = os.path.join(geometry_folder, "temp_geology_input.gpkg")
        temp_raster = os.path.join(geometry_folder, "temp_geology_input.tif")
        temp_definition = os.path.join(
            geometry_folder, "temp_geology_classdefinition.txt"
        )
        dem_crs = self._layer_crs_text(
            QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        )
        definition_backup = self._read_optional_file(classdefinition_path)
        metadata_backup = self._read_optional_file(metadata_path)

        self.log_message("Processing Geology layer through mhm-tools...")
        try:
            lookup_path = local_layer_source(lookup_layer)
            if lookup_path is None:
                lookup_path = materialize_vector_layer(lookup_layer, temp_lookup)

            if isinstance(layer, QgsVectorLayer):
                input_path = local_layer_source(layer)
                if input_path is None:
                    input_path = materialize_vector_layer(layer, temp_vector)
                rasterize_geology_map(
                    input_file=input_path,
                    dem_file=self.filled_dem_path,
                    output_file=output_path,
                    mapping_field=geology_mapping_field,
                    lookup_table=lookup_path,
                    lookup_mapping_field=lookup_mapping_field,
                    input_crs=self._layer_crs_text(layer),
                    dem_crs=dem_crs,
                    log=self.log_message,
                )
                if write_classdefinition:
                    write_geology_classdefinition_file(
                        lookup_table=lookup_path,
                        output_file=temp_definition,
                        log=self.log_message,
                    )
            else:
                input_path = local_layer_source(layer)
                if input_path is None:
                    input_path = self._materialize_geology_raster(layer, temp_raster)
                format_geology_file(
                    input_file=input_path,
                    dem_file=self.filled_dem_path,
                    output_file=output_path,
                    lookup_table=lookup_path,
                    mapping_field=lookup_mapping_field,
                    classdefinition_file=temp_definition,
                    input_crs=self._layer_crs_text(layer),
                    dem_crs=dem_crs,
                    log=self.log_message,
                )

            if not os.path.exists(output_path):
                raise RuntimeError("mhm-tools did not create the geology raster.")
            if write_classdefinition:
                if not os.path.exists(temp_definition):
                    raise RuntimeError(
                        "mhm-tools did not create geology_classdefinition.txt."
                    )
                os.replace(temp_definition, classdefinition_path)
                self.mark_output_prepared(
                    classdefinition_path,
                    name="geology_classdefinition.txt",
                    loaded=False,
                )
                metadata = self._write_geology_metadata_after_processing()
                if not metadata or not os.path.exists(metadata):
                    raise RuntimeError(
                        "pymhm did not create geology_class_metadata.json."
                    )
        except Exception as error:
            self.log_message(f"ERROR preparing geology data: {error}")
            QMessageBox.critical(
                self.dialog,
                "Geology Processing Error",
                f"Could not prepare the geology data:\n{error}",
            )
            self._remove_geology_raster_output(output_path)
            if write_classdefinition:
                self._restore_optional_file(classdefinition_path, definition_backup)
                self._restore_optional_file(metadata_path, metadata_backup)
            self.geology_path = None
            return False
        finally:
            self._cleanup_geology_temporary_inputs(
                temp_lookup, temp_vector, temp_raster, temp_definition
            )

        self.load_layer(output_path, "Geology")
        self.log_message("Geology layer prepared successfully.")
        return True

    def _geology_lookup_input(self):
        """Return the selected lookup layer and mapping field."""
        combo = getattr(self.dialog, "mMapLayerComboBox_geologyLookup", None)
        lookup_layer = combo.currentLayer() if combo is not None else None
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a valid Geology lookup table layer.",
            )
            return None, None

        fields = lookup_layer.fields().names()
        mapping_field = self._selected_lookup_field(
            "comboBox_geologyLookupField", fields
        )
        if not mapping_field:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select the geology lookup field.",
            )
            return None, None
        if not self._required_lookup_field(fields, "GEOLOGY_CLASS"):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Geology lookup table must contain a GEOLOGY_CLASS field.",
            )
            return None, None
        self.log_message(
            f"Using geology lookup fields '{mapping_field}' -> 'GEOLOGY_CLASS'"
        )
        return lookup_layer, mapping_field

    def _materialize_geology_raster(self, layer, output_path):
        """Write a selected QGIS raster to a standalone GeoTIFF."""
        self._remove_geology_raster_output(output_path)
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
            raise RuntimeError("Could not materialize the selected geology raster.")
        return output_path

    def _cleanup_geology_temporary_inputs(
        self, lookup_path, vector_path, raster_path, definition_path
    ):
        """Remove path-based copies made from selected QGIS layers."""
        for vector_output in (lookup_path, vector_path):
            try:
                remove_vector_dataset(vector_output)
            except Exception as error:
                self.log_message(
                    "WARNING: Could not remove temporary geology vector "
                    f"'{vector_output}': {error}"
                )
        try:
            self._remove_geology_raster_output(raster_path)
        except Exception as error:
            self.log_message(
                "WARNING: Could not remove temporary geology raster "
                f"'{raster_path}': {error}"
            )
        try:
            if os.path.exists(definition_path):
                os.remove(definition_path)
        except Exception as error:
            self.log_message(
                "WARNING: Could not remove temporary geology classdefinition "
                f"'{definition_path}': {error}"
            )

    @staticmethod
    def _remove_geology_raster_output(path):
        """Remove a raster and common GDAL sidecars."""
        for candidate in (path, f"{path}.aux.xml", f"{path}.ovr"):
            if os.path.exists(candidate):
                os.remove(candidate)

    def _write_geology_classdefinition_after_processing(self):
        """Regenerate the definition for a cached geology raster."""
        writer = getattr(self, "geology_classification_writer", None)
        if writer is None:
            self.log_message("ERROR: Geology classdefinition writer is not available.")
            return None
        return writer()

    def _write_geology_metadata_after_processing(self):
        """Write plugin-specific geology parameter metadata."""
        writer = getattr(self, "geology_class_metadata_writer", None)
        return writer() if writer is not None else None

    @staticmethod
    def _read_optional_file(path):
        """Read an existing output so a failed refresh can restore it."""
        if not os.path.exists(path):
            return None
        with open(path, "rb") as source:
            return source.read()

    @staticmethod
    def _restore_optional_file(path, content):
        """Restore an old output, or remove a newly created partial output."""
        if content is None:
            if os.path.exists(path):
                os.remove(path)
            return
        with open(path, "wb") as output:
            output.write(content)

    @staticmethod
    def _layer_crs_text(layer):
        """Return a CRS string suitable for assigning missing file metadata."""
        crs = layer.crs() if layer is not None else None
        if crs is None or not crs.isValid():
            return None
        authid = crs.authid()
        return authid or crs.toWkt()

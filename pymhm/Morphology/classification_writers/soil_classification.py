# -*- coding: utf-8 -*-
"""mHM soil classdefinition writer."""
from __future__ import annotations

from ..common import (
    QMessageBox,
    QgsVectorLayer,
    morph_folder,
    os,
    project_geometry_folder,
)
from ..core.base import BaseProcessingMixin
from ..core.soil_sources import (
    local_layer_source,
    materialize_vector_layer,
    remove_vector_dataset,
)
from ...mhm_tools_to_integrate.setup_creation import (
    write_soil_classdefinition_file,
)


class SoilClassificationWriterMixin(BaseProcessingMixin):
    """Write the selected soil lookup table in mHM format."""

    def soil_classdefinition_writer(self, show_success=False) -> str | None:
        """Write ``soil_classdefinition.txt`` through :mod:`mhm_tools`."""
        self.log_message(
            "\n--- Writing Soil Classification Definition File ---")

        if not self.check_prerequisites():
            return None

        lookup_combo = getattr(
            self.dialog, "mMapLayerComboBox_soilLookup", None)
        lookup_layer = (
            lookup_combo.currentLayer() if lookup_combo is not None else None)
        if not lookup_layer:
            self.log_message("ERROR: Soil lookup table layer not selected.")
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a Soil lookup table layer.",
            )
            return None
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Selected soil lookup layer is not valid.")
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select a valid Soil lookup table layer.",
            )
            return None

        output_folder = morph_folder(self.dialog.project_folder)
        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        os.makedirs(output_folder, exist_ok=True)
        os.makedirs(geometry_folder, exist_ok=True)
        output_file = os.path.join(output_folder, "soil_classdefinition.txt")
        temporary_lookup = os.path.join(
            geometry_folder, "temp_soil_classdefinition_lookup.gpkg")
        temporary_output = os.path.join(
            geometry_folder, "temp_soil_classdefinition.txt")

        try:
            lookup_path = local_layer_source(lookup_layer)
            if lookup_path is None:
                lookup_path = materialize_vector_layer(
                    lookup_layer, temporary_lookup)
            result = write_soil_classdefinition_file(
                lookup_table=lookup_path,
                output_file=temporary_output,
                log=self.log_message,
            )
            if not result.is_file():
                raise RuntimeError(
                    "mhm-tools did not create soil_classdefinition.txt.")
            os.replace(temporary_output, output_file)

            self.log_message(
                "Soil classification definition file written successfully: "
                f"{output_file}")
            self.mark_output_prepared(
                output_file,
                name="soil_classdefinition.txt",
                loaded=False,
            )
            if show_success:
                QMessageBox.information(
                    self.dialog,
                    "Success",
                    "Soil classification definition file created successfully."
                    f"\n{output_file}",
                )
            return output_file
        except Exception as error:
            import traceback

            self.log_message(
                "ERROR writing soil classification file: "
                f"{error}\n{traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog,
                "Error",
                f"Error writing soil classification file:\n{error}",
            )
            return None
        finally:
            try:
                remove_vector_dataset(temporary_lookup)
            except Exception as error:
                self.log_message(
                    "WARNING: Could not remove temporary soil lookup "
                    f"'{temporary_lookup}': {error}")
            try:
                if os.path.exists(temporary_output):
                    os.remove(temporary_output)
            except Exception as error:
                self.log_message(
                    "WARNING: Could not remove temporary soil classdefinition "
                    f"'{temporary_output}': {error}")

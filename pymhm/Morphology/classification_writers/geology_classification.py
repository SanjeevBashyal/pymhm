# -*- coding: utf-8 -*-
"""mHM geology classdefinition writer."""
from __future__ import annotations

import json

from ..common import (
    os,
    QMessageBox,
    morph_folder,
    project_geometry_folder,
    QgsVectorLayer,
    NULL,
)
from ..core.base import BaseProcessingMixin


class GeologyClassificationWriterMixin(BaseProcessingMixin):
    """mHM geology classdefinition writer."""

    def geology_classification_writer(self, show_success=False) -> str | None:
        """Write geology_classdefinition.txt from the selected lookup table layer."""
        self.log_message("\n--- Writing Geology Classification Definition File ---")

        if not self.check_prerequisites():
            return None

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
            return None

        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Selected geology lookup layer is not valid.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a valid Geology lookup table layer.")
            return None

        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        output_file = os.path.join(
            morph_output_folder, "geology_classdefinition.txt")
        geometry_output_folder = project_geometry_folder(self.dialog.project_folder)
        os.makedirs(geometry_output_folder, exist_ok=True)
        metadata_file = os.path.join(
            geometry_output_folder, "geology_class_metadata.json")

        try:
            field_names, rows = self._lookup_layer_rows(lookup_layer)
            output_rows = self._geology_classdefinition_rows(field_names, rows)
            if not output_rows:
                raise ValueError("No valid geology classdefinition rows were found.")

            with open(output_file, "w", encoding="utf-8") as output:
                output.write(f"nGeo_Formations  {len(output_rows)}\n")
                output.write("GeoParam(i)   ClassUnit     Karstic      Description\n")
                for row in output_rows:
                    output.write(
                        f"{row['geo_param']:10d}\t"
                        f"{row['class_unit']:10d}     "
                        f"{row['karstic']:10d}      "
                        f"GeoUnit-{row['class_unit']}\n"
                    )
                output.write("!<-END\n")
                output.write("\n")
                output.write("\n")
                output.write("!***********************************\n")
                output.write("! NOTES\n")
                output.write("!***********************************\n")
                output.write("1 = Karstic\n")
                output.write("0 = Non-karstic\n")
                output.write("\n")
                output.write("IMPORTANT ::\n")
                output.write(
                    "   Ordering has to be according to the ordering in mhm_parameter.nml\n")
                output.write("   (namelist: geoparameter)\n")

            self.log_message(
                f"Geology classification definition file written successfully: {output_file}")
            self._write_geology_class_metadata(metadata_file, output_rows)
            self.mark_output_prepared(
                output_file,
                name="geology_classdefinition.txt",
                loaded=False
            )
            if show_success:
                QMessageBox.information(
                    self.dialog, "Success",
                    f"Geology classification definition file created successfully.\n{output_file}")
            return output_file

        except Exception as e:
            import traceback
            self.log_message(
                f"ERROR writing geology classification file: {e}\n{traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing geology classification file:\n{e}")
            return None

    def _lookup_layer_rows(self, lookup_layer):
        """Return lookup field names and feature rows from a QGIS table layer."""
        field_names = lookup_layer.fields().names()
        rows = []
        for feature in lookup_layer.getFeatures():
            rows.append({
                field_name: feature[field_name]
                for field_name in field_names
            })

        self.log_message(
            f"Read {len(rows)} rows from geology lookup layer '{lookup_layer.name()}'.")
        self.log_message(f"Geology lookup fields found: {', '.join(field_names)}")
        return field_names, rows

    def _geology_classdefinition_rows(self, field_names, rows):
        """Build geology classdefinition rows from the template lookup layout."""
        field_lookup = self._geology_field_lookup(field_names)
        geology_class_field = self._required_geology_field(
            field_lookup, "GEOLOGY_CLASS")
        geo_class_field = self._required_geology_field(
            field_lookup, "GEO_CLASS")
        karstic_field = self._required_geology_field(
            field_lookup, "KARSTIC")
        parameter_value_field = self._required_geology_field(
            field_lookup, "PARAMETER_VALUE")

        output_rows = []
        for row_number, row in enumerate(rows, start=2):
            geology_class = self._required_int(
                row.get(geology_class_field), row_number, geology_class_field)
            geo_class = self._required_int(
                row.get(geo_class_field), row_number, geo_class_field)
            karstic = self._required_bool_int(
                row.get(karstic_field), row_number, karstic_field)
            parameter_value = self._required_int(
                row.get(parameter_value_field), row_number, parameter_value_field)

            output_rows.append({
                "geo_param": geo_class,
                "class_unit": geology_class,
                "karstic": karstic,
                "parameter_value": parameter_value,
            })

        return sorted(output_rows, key=lambda item: (
            item["geo_param"], item["class_unit"]))

    def _geology_field_lookup(self, field_names):
        """Map normalized geology lookup field names to original QGIS field names."""
        lookup = {}
        for field_name in field_names:
            lookup.setdefault(
                self._normalise_geology_field_name(field_name),
                field_name,
            )
        return lookup

    def _geology_field(self, field_lookup, field_name):
        """Return a geology lookup field by normalized template name."""
        return field_lookup.get(self._normalise_geology_field_name(field_name))

    def _required_geology_field(self, field_lookup, field_name):
        """Return a required geology lookup field or raise a readable error."""
        field = self._geology_field(field_lookup, field_name)
        if not field:
            raise ValueError(
                f"Geology lookup table is missing required field '{field_name}'.")
        return field

    def _normalise_geology_field_name(self, field_name):
        """Normalize template field names by removing optional marks and units."""
        field_text = str(field_name).strip().lstrip("*").strip()
        if "[" in field_text:
            field_text = field_text.split("[", 1)[0].strip()
        return "".join(char.lower() for char in field_text if char.isalnum())

    def _required_int(self, value, row_number, field_name):
        """Read a required integer value from a lookup row."""
        number = self._required_float(value, row_number, field_name)
        if not number.is_integer():
            raise ValueError(
                f"Row {row_number} has non-integer value '{value}' "
                f"for required field '{field_name}'.")
        return int(number)

    def _required_bool_int(self, value, row_number, field_name):
        """Read a required boolean value and return it as 0 or 1."""
        if self._is_blank(value):
            raise ValueError(
                f"Row {row_number} has an empty value for required field '{field_name}'.")

        if isinstance(value, bool):
            return 1 if value else 0

        value_text = str(value).strip().lower()
        if value_text in ("1", "true", "t", "yes", "y"):
            return 1
        if value_text in ("0", "false", "f", "no", "n"):
            return 0

        raise ValueError(
            f"Row {row_number} has invalid boolean value '{value}' "
            f"for required field '{field_name}'.")

    def _required_float(self, value, row_number, field_name):
        """Read a required numeric value from a lookup row."""
        if self._is_blank(value):
            raise ValueError(
                f"Row {row_number} has an empty value for required field '{field_name}'.")
        try:
            return float(value)
        except (TypeError, ValueError):
            raise ValueError(
                f"Row {row_number} has invalid numeric value '{value}' "
                f"for required field '{field_name}'.")

    def _is_blank(self, value):
        """Return True when a lookup value is empty."""
        return value in (None, NULL, "") or str(value).strip() == ""

    def _write_geology_class_metadata(self, metadata_file, output_rows):
        """Save geology class count and parameter defaults for Configuration."""
        metadata = {
            "version": 1,
            "geology_class_count": len(output_rows),
            "classes": [
                {
                    "geo_param": row["geo_param"],
                    "geology_class": row["class_unit"],
                    "karstic": row["karstic"],
                    "parameter_value": row["parameter_value"],
                }
                for row in output_rows
            ],
        }
        with open(metadata_file, "w", encoding="utf-8") as metadata_output:
            json.dump(metadata, metadata_output, indent=2, sort_keys=True)
        self.mark_output_prepared(
            metadata_file,
            name="geology_class_metadata.json",
            loaded=False,
        )
        self.log_message(
            f"Geology class metadata saved successfully: {metadata_file}")

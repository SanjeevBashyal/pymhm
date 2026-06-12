# -*- coding: utf-8 -*-
"""mHM soil classdefinition writer."""
from __future__ import annotations

from ..common import (
    os,
    QMessageBox,
    morph_folder,
    QgsVectorLayer,
    NULL,
)
from ..core.base import BaseProcessingMixin


class SoilClassificationWriterMixin(BaseProcessingMixin):
    """mHM soil classdefinition writer."""

    def soil_classdefinition_writer(self, show_success=False) -> str | None:
        """Write soil_classdefinition.txt from the selected soil lookup table layer."""
        self.log_message("\n--- Writing Soil Classification Definition File ---")

        if not self.check_prerequisites():
            return None

        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_soilLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )
        if not lookup_layer:
            self.log_message("ERROR: Soil lookup table layer not selected.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table layer.")
            return None

        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Selected soil lookup layer is not valid.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a valid Soil lookup table layer.")
            return None

        morph_output_folder = morph_folder(self.dialog.project_folder)
        os.makedirs(morph_output_folder, exist_ok=True)
        output_file = os.path.join(
            morph_output_folder, "soil_classdefinition.txt")

        try:
            field_names, rows = self._lookup_layer_rows(lookup_layer)
            output_rows = self._soil_classdefinition_rows(field_names, rows)
            if not output_rows:
                raise ValueError("No valid soil horizon rows were found.")

            unique_soil_classes = sorted({
                row["soil_class"] for row in output_rows
            })

            with open(output_file, "w", encoding="utf-8") as output:
                output.write(f"nSoil_Types {len(unique_soil_classes)}\n")
                output.write(
                    "SOIL_NR\tHORIZON\tUD[mm]\tLD[mm]\t"
                    "Clay[%]\tSAND[%]\tBd[gcm-3]\tSilt[%]\n"
                )
                for row in output_rows:
                    output.write(
                        f"{row['soil_class']}\t"
                        f"{row['horizon']}\t"
                        f"{self._format_soil_value(row['upper_depth'])}\t"
                        f"{self._format_soil_value(row['lower_depth'])}\t"
                        f"{self._format_soil_value(row['clay'])}\t"
                        f"{self._format_soil_value(row['sand'])}\t"
                        f"{self._format_soil_value(row['bulk_density'])}\t"
                        f"{self._format_soil_value(row['silt'])}\n"
                    )

            self.log_message(
                f"Soil classification definition file written successfully: {output_file}")
            self.mark_output_prepared(
                output_file,
                name="soil_classdefinition.txt",
                loaded=False
            )
            if show_success:
                QMessageBox.information(
                    self.dialog, "Success",
                    f"Soil classification definition file created successfully.\n{output_file}")
            return output_file

        except Exception as e:
            import traceback
            self.log_message(
                f"ERROR writing soil classification file: {e}\n{traceback.format_exc()}")
            QMessageBox.critical(
                self.dialog, "Error",
                f"Error writing soil classification file:\n{e}")
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
            f"Read {len(rows)} rows from soil lookup layer '{lookup_layer.name()}'.")
        self.log_message(f"Soil lookup fields found: {', '.join(field_names)}")
        return field_names, rows

    def _soil_classdefinition_rows(self, field_names, rows):
        """Build classdefinition rows from either supported soil lookup layout."""
        field_lookup = self._soil_field_lookup(field_names)
        soil_class_field = self._required_soil_field(
            field_lookup, "SOIL_CLASS")

        rowwise_fields = (
            "HORIZON",
            "UPPER_DEPTH",
            "LOWER_DEPTH",
            "CLAY",
            "SAND",
            "SILT",
            "BULK_DENSITY",
        )
        if all(self._soil_field(field_lookup, name) for name in rowwise_fields):
            return self._rowwise_soil_rows(
                rows,
                soil_class_field,
                {
                    name.lower(): self._required_soil_field(field_lookup, name)
                    for name in rowwise_fields
                },
            )

        horizons_field = self._soil_field(field_lookup, "HORIZONS")
        if horizons_field:
            return self._wide_soil_rows(rows, field_lookup, soil_class_field, horizons_field)

        raise ValueError(
            "Soil lookup table must use either the row-per-horizon layout "
            "(HORIZON, UPPER_DEPTH, LOWER_DEPTH, CLAY, SAND, SILT, BULK_DENSITY) "
            "or the wide layout (HORIZONS, DEPTH1, BULK_DENSITY1, CLAY1, SILT1, SAND1, ...)."
        )

    def _rowwise_soil_rows(self, rows, soil_class_field, fields):
        """Parse a lookup table with one row per soil horizon."""
        output_rows = []
        for row_number, row in enumerate(rows, start=2):
            soil_class = self._required_int(
                row.get(soil_class_field), row_number, soil_class_field)
            horizon = self._required_int(
                row.get(fields["horizon"]), row_number, fields["horizon"])
            output_rows.append({
                "soil_class": soil_class,
                "horizon": horizon,
                "upper_depth": self._required_float(
                    row.get(fields["upper_depth"]), row_number, fields["upper_depth"]),
                "lower_depth": self._required_float(
                    row.get(fields["lower_depth"]), row_number, fields["lower_depth"]),
                "clay": self._required_float(
                    row.get(fields["clay"]), row_number, fields["clay"]),
                "sand": self._required_float(
                    row.get(fields["sand"]), row_number, fields["sand"]),
                "silt": self._required_float(
                    row.get(fields["silt"]), row_number, fields["silt"]),
                "bulk_density": self._required_float(
                    row.get(fields["bulk_density"]), row_number, fields["bulk_density"]),
            })

        return sorted(output_rows, key=lambda item: (
            item["soil_class"], item["horizon"]))

    def _wide_soil_rows(self, rows, field_lookup, soil_class_field, horizons_field):
        """Parse a lookup table with one row per soil class and horizon suffixes."""
        output_rows = []
        for row_number, row in enumerate(rows, start=2):
            soil_class = self._required_int(
                row.get(soil_class_field), row_number, soil_class_field)
            horizons = self._required_int(
                row.get(horizons_field), row_number, horizons_field)
            upper_depth = 0.0

            for horizon in range(1, horizons + 1):
                depth_field = self._required_soil_field(
                    field_lookup, f"DEPTH{horizon}")
                bulk_density_field = self._required_soil_field(
                    field_lookup, f"BULK_DENSITY{horizon}")
                clay_field = self._required_soil_field(
                    field_lookup, f"CLAY{horizon}")
                silt_field = self._required_soil_field(
                    field_lookup, f"SILT{horizon}")
                sand_field = self._required_soil_field(
                    field_lookup, f"SAND{horizon}")

                lower_depth = self._required_float(
                    row.get(depth_field), row_number, depth_field)
                output_rows.append({
                    "soil_class": soil_class,
                    "horizon": horizon,
                    "upper_depth": upper_depth,
                    "lower_depth": lower_depth,
                    "clay": self._required_float(
                        row.get(clay_field), row_number, clay_field),
                    "sand": self._required_float(
                        row.get(sand_field), row_number, sand_field),
                    "silt": self._required_float(
                        row.get(silt_field), row_number, silt_field),
                    "bulk_density": self._required_float(
                        row.get(bulk_density_field), row_number, bulk_density_field),
                })
                upper_depth = lower_depth

        return sorted(output_rows, key=lambda item: (
            item["soil_class"], item["horizon"]))

    def _soil_field_lookup(self, field_names):
        """Map normalized soil lookup field names to original QGIS field names."""
        lookup = {}
        for field_name in field_names:
            lookup.setdefault(
                self._normalise_soil_field_name(field_name),
                field_name,
            )
        return lookup

    def _soil_field(self, field_lookup, field_name):
        """Return a soil lookup field by normalized template name."""
        return field_lookup.get(self._normalise_soil_field_name(field_name))

    def _required_soil_field(self, field_lookup, field_name):
        """Return a required soil lookup field or raise a readable error."""
        field = self._soil_field(field_lookup, field_name)
        if not field:
            raise ValueError(f"Soil lookup table is missing required field '{field_name}'.")
        return field

    def _normalise_soil_field_name(self, field_name):
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

    def _format_soil_value(self, value):
        """Format numbers like the soil classdefinition template."""
        if value is None:
            return ""
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)

        if number.is_integer():
            return str(int(number))
        return f"{number:.6g}"

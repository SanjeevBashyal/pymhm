# -*- coding: utf-8 -*-
"""Soil lookup table parsing."""
from ..common import (
    QMessageBox,
    QgsVectorLayer,
)
from .lookup_fields import LookupFieldMixin


class SoilLookupMixin(LookupFieldMixin):
    """Soil lookup table parsing."""

    def load_soil_lookup_table(self):
        """Load soil lookup table from the selected QGIS table layer."""
        self.soil_lookup_key_field = None
        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_soilLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )
        if not lookup_layer:
            self.log_message(
                "ERROR: Soil lookup table layer not selected.")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select a Soil lookup table layer.")
            return None

        return self._read_soil_lookup_layer(lookup_layer)

    def _read_soil_lookup_layer(self, lookup_layer):
        """Read soil lookup mappings from a selected QGIS table layer."""
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Soil lookup layer is not valid.")
            return None

        field_names = lookup_layer.fields().names()
        key_field = self._selected_lookup_field(
            "comboBox_soilLookupField", field_names)
        if not key_field:
            self.log_message(
                "ERROR: Please select a soil lookup field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Please select the soil lookup field that maps the soil layer to the lookup table.")
            return None

        class_field = self._required_lookup_field(field_names, "CLASS")
        if not class_field:
            self.log_message("ERROR: Soil lookup layer must contain a CLASS field.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            QMessageBox.warning(
                self.dialog, "Input Error",
                "Soil lookup layer must contain a CLASS field.")
            return None

        self.soil_lookup_key_field = key_field
        self.log_message(
            f"Using soil lookup fields '{key_field}' -> '{class_field}'")

        lookup_mapping = {}
        for feature_index, feature in enumerate(lookup_layer.getFeatures(), start=1):
            try:
                lookup_key = self._normalise_lookup_key(feature[key_field])
                if not lookup_key:
                    self.log_message(
                        f"WARNING: Soil lookup row {feature_index} has an empty key. Skipping.")
                    continue

                class_code = self._coerce_lookup_int(feature[class_field])
                if class_code is None or class_code <= 0:
                    self.log_message(
                        f"WARNING: Soil lookup row {feature_index} has invalid CLASS '{feature[class_field]}'. Skipping.")
                    continue

                lookup_mapping[lookup_key] = class_code
                if len(lookup_mapping) <= 10:
                    self.log_message(f"  Mapped: {lookup_key} -> {class_code}")
            except Exception as e:
                self.log_message(
                    f"WARNING: Error parsing soil lookup feature {feature_index}: {e}")

        if lookup_mapping:
            self.log_message(
                f"Loaded {len(lookup_mapping)} entries from soil lookup layer.")
            return lookup_mapping

        self.log_message("ERROR: Soil lookup layer is empty or could not be parsed.")
        QMessageBox.warning(
            self.dialog, "Input Error",
            "Soil lookup layer is empty or could not be parsed.")
        return None

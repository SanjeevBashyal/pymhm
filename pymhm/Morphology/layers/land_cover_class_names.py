# -*- coding: utf-8 -*-
"""Raw land-cover class-name lookup for reporting outputs."""
from __future__ import annotations

from typing import Any

from ..common import (
    os,
    csv,
    QgsVectorLayer,
    NULL,
)
from .lookup_fields import LookupFieldMixin


class LandCoverClassNameMixin(LookupFieldMixin):
    """Raw land-cover class-name lookup for reporting outputs."""

    def _read_land_cover_class_names(self) -> dict[int, str]:
        """Read raw land-cover class value to class-name mappings."""
        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_landCoverLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )

        if lookup_layer:
            class_names = self._read_land_cover_class_names_from_layer(
                lookup_layer)
            if class_names:
                return class_names

        lookup_file_widget = getattr(
            self.dialog, "lineEdit_land_cover_lookup", None)
        lookup_file = lookup_file_widget.text() if lookup_file_widget else ""
        if lookup_file and os.path.exists(lookup_file):
            return self._read_land_cover_class_names_from_csv(lookup_file)

        return {}

    def _read_land_cover_class_names_from_layer(
            self,
            lookup_layer: Any) -> dict[int, str]:
        """Read class names from the selected tabular lookup layer."""
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            return {}

        field_names = lookup_layer.fields().names()
        normalised_fields = {
            self._normalise_lookup_field_name(field_name): field_name
            for field_name in field_names
        }

        grid_field = self._selected_lookup_field(
            "comboBox_landCoverLookupField", field_names)
        if not grid_field:
            grid_field = self._first_lookup_field(
                normalised_fields,
                ("gridvalue", "grid", "value", "id", "classvalue", "landcoverid")
            )
        name_field = self._first_lookup_field(
            normalised_fields,
            (
                "classname",
                "name",
                "label",
                "description",
                "landcoverclass",
                "landcovername",
                "landcover",
                "category",
                "type",
                "class"
            )
        )

        if not grid_field or not name_field:
            self.log_message(
                "WARNING: Land cover lookup layer has no class-name fields for band details.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            return {}

        class_names = {}
        for feature in lookup_layer.getFeatures():
            try:
                class_value = self._coerce_lookup_int(feature[grid_field])
                raw_class_name = feature[name_field]
                if raw_class_name in (None, NULL, ""):
                    continue
                class_name = str(raw_class_name).strip()
                if class_value is None or not class_name:
                    continue
                class_names[class_value] = class_name
            except Exception as e:
                self.log_message(
                    f"WARNING: Error parsing land cover class-name feature: {e}")

        return class_names

    def _read_land_cover_class_names_from_csv(
            self,
            lookup_file: str) -> dict[int, str]:
        """Read class names from the older CSV lookup file."""
        try:
            with open(lookup_file, "r", newline="", encoding="utf-8") as csv_file:
                reader = csv.DictReader(csv_file)
                field_names = reader.fieldnames or []
                normalised_fields = {
                    self._normalise_lookup_field_name(field_name): field_name
                    for field_name in field_names
                }

                grid_field = self._first_lookup_field(
                    normalised_fields,
                    ("gridvalue", "grid", "value", "id", "classvalue", "landcoverid")
                )
                name_field = self._first_lookup_field(
                    normalised_fields,
                    (
                        "classname",
                        "name",
                        "label",
                        "description",
                        "landcoverclass",
                        "landcovername",
                        "landcover",
                        "category",
                        "type",
                        "class"
                    )
                )

                if not grid_field or not name_field:
                    return {}

                class_names = {}
                for row in reader:
                    try:
                        class_value = self._coerce_lookup_int(row.get(grid_field))
                        raw_class_name = row.get(name_field, "")
                        if raw_class_name in (None, ""):
                            continue
                        class_name = str(raw_class_name).strip()
                        if class_value is None or not class_name:
                            continue
                        class_names[class_value] = class_name
                    except Exception as e:
                        self.log_message(
                            f"WARNING: Error parsing land cover class-name row: {e}")

                return class_names
        except Exception as e:
            self.log_message(
                f"WARNING: Could not read land cover class names from CSV: {e}")
            return {}

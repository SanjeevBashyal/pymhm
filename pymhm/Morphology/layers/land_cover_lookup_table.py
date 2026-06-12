# -*- coding: utf-8 -*-
"""Land-cover lookup table loading for mHM class reclassification."""
from __future__ import annotations

from typing import Any

from ..common import (
    os,
    csv,
    QgsVectorLayer,
)
from .lookup_fields import LookupFieldMixin


class LandCoverLookupTableMixin(LookupFieldMixin):
    """Land-cover lookup table loading for mHM class reclassification."""

    def load_land_cover_lookup_table(self) -> dict[int, int]:
        """
        Load land cover lookup table from a selected tabular layer, an older
        file-path widget, or use the default mapping.

        Returns:
            dict: Mapping of grid values to type_int (1=Forest, 2=Impervious, 3=Pervious)
        """
        # Default lookup table
        default_lookup = {
            1: 3,   # Waterbody -> Pervious
            2: 2,   # Glacier -> Impervious
            3: 3,   # Snow -> Pervious
            4: 1,   # Forest -> Forest
            5: 3,   # Riverbed -> Pervious
            6: 2,   # Built-up area -> Impervious
            7: 3,   # Cropland -> Pervious
            8: 3,   # Bare soil -> Pervious
            9: 2,   # Bare rock -> Impervious
            10: 3,  # Grassland -> Pervious
            11: 3   # Other wooded land -> Pervious
        }

        lookup_layer_combo = getattr(
            self.dialog, "mMapLayerComboBox_landCoverLookup", None)
        lookup_layer = (
            lookup_layer_combo.currentLayer()
            if lookup_layer_combo is not None else None
        )

        if lookup_layer:
            lookup_mapping = self._read_land_cover_lookup_layer(lookup_layer)
            if lookup_mapping:
                self.log_message(
                    f"Loaded {len(lookup_mapping)} entries from land cover lookup layer.")
                return lookup_mapping

            self.log_message(
                "Land cover lookup layer could not be parsed. Using default mapping.")
            return default_lookup

        lookup_file_widget = getattr(
            self.dialog, "lineEdit_land_cover_lookup", None)
        lookup_file = lookup_file_widget.text() if lookup_file_widget else ""

        if not lookup_file or not os.path.exists(lookup_file):
            self.log_message(
                "No land cover lookup layer/file provided. Using default mapping.")
            return default_lookup

        # Try to read lookup table from file (CSV format expected)
        try:
            import csv
            lookup_mapping = {}
            with open(lookup_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        grid_value = int(row.get('Grid value', 0))
                        type_str = row.get('Type', '').strip()

                        # Map type string to type_int
                        if type_str == 'Forest':
                            type_int = 1
                        elif type_str == 'Impervious':
                            type_int = 2
                        elif type_str == 'Pervious':
                            type_int = 3
                        else:
                            self.log_message(
                                f"WARNING: Unknown type '{type_str}' in lookup table. Skipping.")
                            continue

                        lookup_mapping[grid_value] = type_int
                    except (ValueError, KeyError) as e:
                        self.log_message(
                            f"WARNING: Error parsing lookup table row: {e}")
                        continue

            if lookup_mapping:
                self.log_message(
                    f"Loaded {len(lookup_mapping)} entries from lookup table.")
                return lookup_mapping
            else:
                self.log_message(
                    "Lookup table is empty. Using default mapping.")
                return default_lookup
        except Exception as e:
            self.log_message(
                f"ERROR reading lookup table: {e}. Using default mapping.")
            return default_lookup

    def _read_land_cover_lookup_layer(
            self,
            lookup_layer: Any) -> dict[int, int] | None:
        """Read land-cover lookup mappings from a tabular/vector layer."""
        if not isinstance(lookup_layer, QgsVectorLayer) or not lookup_layer.isValid():
            self.log_message("ERROR: Land cover lookup layer is not valid.")
            return None

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
        type_field = self._first_lookup_field(
            normalised_fields,
            ("type", "typeint", "mhmtype", "mhmclass", "class", "landcovertype")
        )

        if not grid_field or not type_field:
            self.log_message(
                "ERROR: Land cover lookup layer must contain grid value and type fields.")
            self.log_message(f"Fields found: {', '.join(field_names)}")
            return None

        lookup_mapping = {}
        for feature in lookup_layer.getFeatures():
            try:
                grid_value = self._coerce_lookup_int(feature[grid_field])
                type_int = self._coerce_land_cover_type(feature[type_field])
                if grid_value is None or type_int is None:
                    continue
                lookup_mapping[grid_value] = type_int
            except Exception as e:
                self.log_message(
                    f"WARNING: Error parsing land cover lookup feature: {e}")

        return lookup_mapping or None

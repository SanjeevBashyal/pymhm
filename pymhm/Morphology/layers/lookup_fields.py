# -*- coding: utf-8 -*-
"""Shared lookup-table field matching and value coercion helpers."""
from ..common import NULL


class LookupFieldMixin:
    """Shared lookup-table field matching and value coercion helpers."""

    def _normalise_lookup_field_name(self, field_name):
        """Normalise a lookup-table field name for case-insensitive matching."""
        field_text = str(field_name).strip().lstrip("*").strip()
        if "[" in field_text:
            field_text = field_text.split("[", 1)[0].strip()
        return "".join(
            char.lower() for char in field_text if char.isalnum())

    def _first_lookup_field(self, normalised_fields, candidates):
        """Return the original field name for the first matching candidate."""
        for candidate in candidates:
            if candidate in normalised_fields:
                return normalised_fields[candidate]
        return None

    def _selected_lookup_field(self, combo_name, field_names):
        """Return a selected lookup-table field if it exists in field_names."""
        combo_box = getattr(self.dialog, combo_name, None)
        if combo_box is None:
            return None

        selected_field = combo_box.currentText().strip()
        if not selected_field:
            return None

        if selected_field in field_names:
            return selected_field

        normalised_fields = {
            self._normalise_lookup_field_name(field_name): field_name
            for field_name in field_names
        }
        return normalised_fields.get(
            self._normalise_lookup_field_name(selected_field))

    def _required_lookup_field(self, field_names, required_name):
        """Return required_name from field_names using normalised exact matching."""
        normalised_fields = {
            self._normalise_lookup_field_name(field_name): field_name
            for field_name in field_names
        }
        return normalised_fields.get(
            self._normalise_lookup_field_name(required_name))

    def _normalise_lookup_key(self, value):
        """Convert lookup keys to stable strings for table-to-layer matching."""
        if value in (None, NULL, ""):
            return ""

        value_text = str(value).strip()
        if not value_text:
            return ""

        try:
            value_number = float(value_text)
        except ValueError:
            return value_text

        if value_number.is_integer():
            return str(int(value_number))
        return value_text

    def _coerce_lookup_int(self, value):
        """Convert a table value to int, treating NULL/blank as missing."""
        if value in (None, NULL, ""):
            return None
        number = float(value)
        if not number.is_integer():
            return None
        return int(number)

    def _coerce_land_cover_type(self, value):
        """Convert land-cover type strings or numeric codes to mHM class IDs."""
        if value in (None, NULL, ""):
            return None

        type_text = str(value).strip()
        type_key = type_text.lower()
        type_lookup = {
            "forest": 1,
            "impervious": 2,
            "pervious": 3,
        }
        if type_key in type_lookup:
            return type_lookup[type_key]

        try:
            return int(float(type_text))
        except ValueError:
            self.log_message(
                f"WARNING: Unknown land cover type '{type_text}'. Skipping.")
            return None

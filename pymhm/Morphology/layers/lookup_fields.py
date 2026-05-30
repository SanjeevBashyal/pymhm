# -*- coding: utf-8 -*-
"""Shared lookup-table field matching and value coercion helpers."""
from ..common import NULL


class LookupFieldMixin:
    """Shared lookup-table field matching and value coercion helpers."""

    def _normalise_lookup_field_name(self, field_name):
        """Normalise a lookup-table field name for case-insensitive matching."""
        return "".join(
            char.lower() for char in str(field_name) if char.isalnum())

    def _first_lookup_field(self, normalised_fields, candidates):
        """Return the original field name for the first matching candidate."""
        for candidate in candidates:
            if candidate in normalised_fields:
                return normalised_fields[candidate]
        return None

    def _coerce_lookup_int(self, value):
        """Convert a table value to int, treating NULL/blank as missing."""
        if value in (None, NULL, ""):
            return None
        return int(float(value))

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

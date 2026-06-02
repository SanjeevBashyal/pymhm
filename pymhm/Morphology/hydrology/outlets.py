# -*- coding: utf-8 -*-
"""Gauge outlet helpers for the Hydrology page."""
from __future__ import annotations

STATION_ID_FIELD = "STATION_ID"


class StationIdError(ValueError):
    """Raised when a pour-point layer does not provide usable station IDs."""


def find_station_id_field(layer) -> str:
    """Return the station id field name, accepting case-insensitive matches."""
    if not layer or not layer.isValid():
        raise StationIdError("Please select a valid pour points layer.")

    field_names = layer.fields().names()
    if STATION_ID_FIELD in field_names:
        return STATION_ID_FIELD

    for field_name in field_names:
        if field_name.upper() == STATION_ID_FIELD:
            return field_name

    raise StationIdError(
        f"The pour points layer must contain a {STATION_ID_FIELD} column.")


def station_id_text(value) -> str:
    """Convert a station id attribute to a stable filename-safe text value."""
    if value is None:
        return ""

    text = str(value).strip()
    if not text or text.upper() == "NULL":
        return ""

    try:
        numeric = float(text)
        if numeric.is_integer():
            return str(int(numeric))
    except Exception:
        pass

    return text


def station_id_int(value) -> int:
    """Convert a station id attribute to the integer value burned in idgauges."""
    text = station_id_text(value)
    if not text:
        raise StationIdError("Empty STATION_ID value found.")

    try:
        return int(float(text))
    except Exception as e:
        raise StationIdError(
            f"STATION_ID value '{text}' cannot be burned into a raster as an integer."
        ) from e


def station_ids_from_layer(layer, unique: bool = True) -> list[str]:
    """Return station IDs from a pour-point layer in feature order."""
    field_name = find_station_id_field(layer)
    station_ids = []
    seen = set()

    for feature in layer.getFeatures():
        station_id = station_id_text(feature.attribute(field_name))
        if not station_id:
            raise StationIdError("A pour point feature has an empty STATION_ID.")
        if unique:
            if station_id in seen:
                raise StationIdError(
                    f"Duplicate STATION_ID '{station_id}' found in pour points.")
            seen.add(station_id)
        station_ids.append(station_id)

    return station_ids


def feature_count_text(layer) -> str:
    """Return a display value for the Hydrology page gauged outlet count."""
    if not layer or not layer.isValid():
        return "0"
    try:
        return str(int(layer.featureCount()))
    except Exception:
        return "0"


class OutletCountMixin:
    """Hydrology page gauged outlet count display."""

    def update_gauged_outlet_count(self, layer=None):
        """Show the selected pour-point feature count in the Hydrology page."""
        if layer is None:
            layer = self.dialog.mMapLayerComboBox_pour_points.currentLayer()

        count = feature_count_text(layer)
        if hasattr(self.dialog, "label_numberOfGaugedOutletsValue"):
            self.dialog.label_numberOfGaugedOutletsValue.setText(count)
        return count

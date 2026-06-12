# -*- coding: utf-8 -*-
"""Domain and class discovery for Configuration dialogs."""

from __future__ import annotations

import json
import os
import re
from typing import Any

from ..project_layout import geometry_folder, morph_folder


def _field_name(layer: Any, wanted: str) -> str | None:
    """Return a field name by case-insensitive lookup."""
    if not layer:
        return None
    wanted = wanted.lower()
    for field in layer.fields():
        if field.name().lower() == wanted:
            return field.name()
    return None


def _truthy(value: Any) -> bool:
    """Return True for common boolean-ish attribute values."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in ("1", "true", "t", "yes", "y")


def domain_infos(dialog: Any) -> list[dict[str, int | str]]:
    """
    Return domain labels from pour points where IS_DOMAIN is true.

    Falls back to one domain labelled "1" when no valid domain feature can be
    discovered.
    """
    layer = None
    if hasattr(dialog, "mMapLayerComboBox_pour_points"):
        layer = dialog.mMapLayerComboBox_pour_points.currentLayer()
    if not layer or not layer.isValid():
        return [{"station_id": "1", "label": "1", "index": 1}]

    station_field = _field_name(layer, "STATION_ID")
    is_domain_field = _field_name(layer, "IS_DOMAIN")
    if not station_field or not is_domain_field:
        return [{"station_id": "1", "label": "1", "index": 1}]

    domains = []
    for feature in layer.getFeatures():
        if not _truthy(feature.attribute(is_domain_field)):
            continue
        station_id = str(feature.attribute(station_field)).strip()
        if not station_id:
            continue
        domains.append({
            "station_id": station_id,
            "label": station_id,
            "index": len(domains) + 1,
        })

    return domains or [{"station_id": "1", "label": "1", "index": 1}]


def domain_count(dialog: Any) -> int:
    """Return the current domain count."""
    return len(domain_infos(dialog))


def geology_class_count(dialog: Any, default: int = 16) -> int:
    """Return the number of geological classes from prepared morphology data."""
    return len(geology_class_rows(dialog, default))


def geology_class_rows(dialog: Any, default: int = 16) -> list[dict[str, Any]]:
    """Return GeoParam-indexed geology class rows for Configuration."""
    metadata = geology_class_metadata(dialog)
    rows = []
    for row in metadata.get("classes", []):
        if not isinstance(row, dict):
            continue
        try:
            geo_param = int(row.get("geo_param"))
            geology_class = int(row.get("geology_class"))
        except (TypeError, ValueError):
            continue
        rows.append({
            "geo_param": geo_param,
            "geology_class": geology_class,
            "karstic": row.get("karstic"),
            "parameter_value": row.get("parameter_value"),
        })

    if rows:
        return sorted(rows, key=lambda item: (
            item["geo_param"], item["geology_class"]))

    count = _geology_classdefinition_count(dialog, default)
    return [
        {
            "geo_param": index,
            "geology_class": index,
            "karstic": None,
            "parameter_value": None,
        }
        for index in range(1, count + 1)
    ]


def geology_class_metadata(dialog: Any) -> dict[str, Any]:
    """Return saved geology class metadata for Configuration defaults."""
    if not getattr(dialog, "project_folder", None):
        return {}

    metadata_path = os.path.join(
        geometry_folder(dialog.project_folder),
        "geology_class_metadata.json",
    )
    if not os.path.exists(metadata_path):
        return {}

    try:
        with open(metadata_path, "r", encoding="utf-8") as metadata_file:
            metadata = json.load(metadata_file)
        if isinstance(metadata, dict):
            classes = metadata.get("classes")
            if isinstance(classes, list):
                return metadata
    except Exception:
        pass
    return {}


def geology_parameter_values(dialog: Any) -> dict[str, Any]:
    """Return GeoParam-indexed PARAMETER_VALUE defaults from geology metadata."""
    values = {}
    for row in geology_class_rows(dialog):
        geo_param = row.get("geo_param")
        parameter_value = row.get("parameter_value")
        try:
            geo_param = int(geo_param)
        except (TypeError, ValueError):
            continue
        if parameter_value is None:
            continue
        values[str(geo_param)] = parameter_value
    return values


def _geology_classdefinition_count(dialog: Any, default: int = 16) -> int:
    """Return geology class count from the classdefinition text file."""
    if not getattr(dialog, "project_folder", None):
        return default

    classdefinition = os.path.join(
        morph_folder(dialog.project_folder),
        "geology_classdefinition.txt",
    )
    if os.path.exists(classdefinition):
        try:
            with open(classdefinition, "r", encoding="utf-8") as file_obj:
                for line in file_obj:
                    match = re.search(r"nGeo_Formations\s+(\d+)", line)
                    if match:
                        return max(1, int(match.group(1)))
        except Exception:
            pass

    return default

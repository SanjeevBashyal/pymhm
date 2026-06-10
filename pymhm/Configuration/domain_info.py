# -*- coding: utf-8 -*-
"""Domain and class discovery for Configuration dialogs."""

from __future__ import annotations

import os
import re
from typing import Any

from ..project_layout import morph_folder


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

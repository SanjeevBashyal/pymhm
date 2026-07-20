# -*- coding: utf-8 -*-
"""Convert dialog values into converter-ready namelist profile values."""
from __future__ import annotations

from typing import Any

from ..project_layout import version_key
from .namelist import canonical_name
from .version_compat import render_values_for_version


ConfigValues = dict[str, dict[str, Any]]


def namelist_values(
        version_text: str,
        kind: str,
        editor_values: ConfigValues,
        pages: list[dict[str, Any]],
        dialog: Any | None = None) -> ConfigValues:
    """Return exact names and JSON array shapes for one file profile."""
    if kind == "mhm" and version_key(version_text) == "v5.13":
        return render_values_for_version(
            version_text, kind, editor_values, dialog)

    values = schema_values(editor_values, pages, version_text, kind)
    if version_key(version_text) != "v5.13":
        return values
    if kind == "parameters":
        return _v513_parameter_values(values)
    if kind == "outputs":
        return _v513_output_values(values)
    return values


def schema_values(
        editor_values: ConfigValues,
        pages: list[dict[str, Any]],
        version_text: str,
        kind: str) -> ConfigValues:
    """Restore exact schema names and reshape internal indexed rows."""
    result: ConfigValues = {}
    for page in pages:
        block = page["block"]
        source = _block(editor_values, block)
        block_values: dict[str, Any] = {}
        for name, prop_schema in page["schema"].get("properties", {}).items():
            value = _property_value(
                source, name, prop_schema, version_text, kind)
            if value is not None and value != []:
                block_values[name] = value
        result[block] = block_values
    return result


def _block(values: ConfigValues, name: str) -> dict[str, Any]:
    key = canonical_name(name)
    for candidate, block_values in (values or {}).items():
        if canonical_name(candidate) == key and isinstance(block_values, dict):
            return block_values
    return {}


def _property_value(
        values: dict[str, Any],
        name: str,
        prop_schema: dict[str, Any],
        version_text: str,
        kind: str) -> Any:
    key = canonical_name(name)
    for candidate, value in values.items():
        if canonical_name(candidate) == key:
            return value

    suffix = "class" if key == "geoparam" else "domain"
    prefix = f"{key}__{suffix}"
    indexed = []
    for candidate, value in values.items():
        if not candidate.startswith(prefix):
            continue
        index = candidate[len(prefix):]
        if index.isdigit():
            indexed.append((int(index), value))
    if not indexed:
        return None

    rows = [value for _, value in sorted(indexed)]
    if key == "geoparam" and version_key(version_text) == "v5.13":
        return rows
    shape = prop_schema.get("x-fortran-shape")
    if isinstance(shape, (list, tuple)) and len(shape) >= 2:
        return _transpose(rows)
    return rows


def _transpose(rows: list[Any]) -> list[list[Any]]:
    if not rows or not all(isinstance(row, list) for row in rows):
        return rows
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("indexed namelist rows must have equal lengths")
    return [[row[column] for row in rows] for column in range(width)]


def _v513_parameter_values(values: ConfigValues) -> ConfigValues:
    aliases = {
        "directrunoff1": "directRunoff1",
        "petm1": "PETminus1",
        "petm2": "PET0",
    }
    return {
        aliases.get(canonical_name(block), block): block_values
        for block, block_values in values.items()
    }


def _value(values: dict[str, Any], name: str, default: Any = None) -> Any:
    key = canonical_name(name)
    for candidate, value in values.items():
        if canonical_name(candidate) == key:
            return value
    return default


def _v513_output_values(values: ConfigValues) -> ConfigValues:
    source = _block(values, "output_mhm")
    fields = {
        1: "out_interception",
        2: "out_snowpack",
        3: "out_SWC",
        4: "out_SM",
        5: "out_SM_all",
        6: "out_sealedSTW",
        7: "out_unsatSTW",
        8: "out_satSTW",
        9: "out_PET",
        10: "out_aET_all",
        11: "out_Q",
        12: "out_QD",
        13: "out_QIf",
        14: "out_QIs",
        15: "out_QB",
        16: "out_recharge",
        17: "out_soil_infil",
        18: "out_neutrons",
        19: "out_aET_layer",
        20: "out_preEffect",
        21: "out_Qsm",
    }
    return {
        "NLoutputResults": {
            "output_deflate_level": _value(source, "output_deflate_level", 6),
            "output_double_precision": _value(
                source, "output_double_precision", False),
            "output_time_reference": _value(
                source, "output_time_reference", 0),
            "timeStep_model_outputs": _value(source, "output_frequency", -2),
            "outputFlxState": [
                bool(_value(source, fields[index], False))
                for index in range(1, 22)
            ],
        }
    }

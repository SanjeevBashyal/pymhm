"""Plugin-specific geology parameter metadata."""

from __future__ import annotations

import json
import os
from pathlib import Path


def _field(columns, requested):
    normalized = _normalize(requested)
    matches = [
        column
        for column in columns
        if not str(column).strip().startswith("*")
        and _normalize(column) == normalized
    ]
    if len(matches) != 1:
        available = ", ".join(str(column) for column in columns)
        raise ValueError(
            f"Geology lookup field {requested!r} was not found uniquely. "
            f"Available fields: {available or '<none>'}."
        )
    return matches[0]


def _optional_field(columns, *names):
    """Resolve the first available non-starred optional geology field."""
    for name in names:
        try:
            return _field(columns, name)
        except ValueError:
            pass
    return None


def _normalize(value):
    text = str(value).strip().lstrip("*").split("[", 1)[0]
    return "".join(char.lower() for char in text if char.isalnum())


def _integer(value, field, row):
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Row {row} has invalid value {value!r} for {field!r}."
        ) from error
    if not number.is_integer():
        raise ValueError(f"Row {row} has non-integer value {value!r} for {field!r}.")
    return int(number)


def _boolean(value, field, row):
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return 1
    if text in {"0", "false", "f", "no", "n"}:
        return 0
    raise ValueError(f"Row {row} has invalid boolean value {value!r} for {field!r}.")


def write_geology_metadata(lookup_table, class_field, output_file):
    """Write metadata consumed by mhm_qgis's geology parameter configuration."""
    from .applications.mhm_tools_handler import read_categorical_lookup_table

    table = read_categorical_lookup_table(lookup_table)
    class_column = _field(table.columns, class_field)
    geo_column = _optional_field(table.columns, "GEO_CLASS", "GEO_ID")
    geo_column = geo_column or class_column
    karst_column = _field(table.columns, "KARSTIC")
    parameter_column = _field(table.columns, "PARAMETER_VALUE")
    rows = []
    for row_number, (_, row) in enumerate(table.iterrows(), start=2):
        rows.append(
            {
                "geo_param": _integer(row[geo_column], geo_column, row_number),
                "geology_class": _integer(
                    row[class_column], class_column, row_number
                ),
                "karstic": _boolean(row[karst_column], karst_column, row_number),
                "parameter_value": _integer(
                    row[parameter_column], parameter_column, row_number
                ),
            }
        )
    rows.sort(key=lambda row: (row["geo_param"], row["geology_class"]))
    metadata = {
        "version": 1,
        "geology_class_count": len(rows),
        "classes": rows,
    }
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(f"{output.suffix}.tmp")
    temporary.write_text(
        json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    os.replace(temporary, output)
    return output


__all__ = ["write_geology_metadata"]

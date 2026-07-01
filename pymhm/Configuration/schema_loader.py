# -*- coding: utf-8 -*-
"""Load YAML namelist schemas."""
from __future__ import annotations

import glob
import os
from typing import Any

from .constants import CONFIG_ORDER, OUTPUT_ORDER
from .paths import schemas_dir


def load_yaml_package() -> Any:
    """Import PyYAML lazily so the plugin can report a useful error."""
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError(
            "PyYAML is required to read the namelist schema files. "
            "Install the 'yaml' Python package in the QGIS Python "
            "environment and try again."
        ) from exc
    return yaml


def load_schema(path: str) -> dict[str, Any]:
    """Load one YAML schema using yaml.safe_load."""
    yaml = load_yaml_package()
    with open(path, "r", encoding="utf-8") as schema_file:
        data = yaml.safe_load(schema_file) or {}
    data["_schema_path"] = path
    return data


def load_schema_group(
        prefix: str,
        version_text: str | None = None) -> dict[str, dict[str, Any]]:
    """Load all schemas with a prefix and index them by fortran namelist."""
    pattern = os.path.join(schemas_dir(version_text), f"{prefix}_*.yml")
    schemas = {}
    for path in sorted(glob.glob(pattern)):
        schema = load_schema(path)
        namelist = schema.get("x-fortran-namelist")
        if namelist:
            schemas[namelist.lower()] = schema
    return schemas


def ordered_schemas(
        prefix: str,
        order: list[str] | tuple[str, ...] | None = None,
        version_text: str | None = None) -> list[dict[str, Any]]:
    """Load schemas and return them in a requested namelist order."""
    schemas = load_schema_group(prefix, version_text)
    ordered = []
    if order:
        for name in order:
            schema = schemas.pop(name.lower(), None)
            if schema:
                ordered.append(schema)
    ordered.extend(schemas[name] for name in sorted(schemas))
    return ordered


def config_schemas(version_text: str | None = None) -> list[dict[str, Any]]:
    """Return mHM configuration schemas in the requested page order."""
    return ordered_schemas("config", CONFIG_ORDER, version_text)


def output_schemas(version_text: str | None = None) -> list[dict[str, Any]]:
    """Return output schemas in template order."""
    return ordered_schemas("output", OUTPUT_ORDER, version_text)


def parameter_schema_lookup(version_text: str | None = None) -> dict[str, dict[str, Any]]:
    """Return parameter schemas indexed by namelist name."""
    return load_schema_group("parameter", version_text)

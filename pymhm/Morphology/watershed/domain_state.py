# -*- coding: utf-8 -*-
"""QGIS-free persistence helpers for watershed domain delineation."""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Mapping, Union

from ...project_layout import workspace_folder


STATE_FILENAME = "pymhm_domain_delineation_state.json"
STATE_VERSION = 1
DEM_DOMAIN_ID = "__dem__"
PathInput = Union[str, Path]


def default_state() -> Dict[str, Any]:
    """Return a new empty domain-delineation state."""
    return {
        "version": STATE_VERSION,
        "pour_points_source": "",
        "outlet_id_field": "",
        "dem_domain": False,
        "outlets": {},
    }


def state_path(project_folder: PathInput) -> Path:
    """Return the state path inside the plugin-owned workspace."""
    return Path(workspace_folder(project_folder)) / STATE_FILENAME


def _normalized_state(value: object) -> Dict[str, Any]:
    state = default_state()
    if not isinstance(value, Mapping):
        return state

    for key in ("pour_points_source", "outlet_id_field"):
        item = value.get(key)
        if isinstance(item, str):
            state[key] = item

    if isinstance(value.get("dem_domain"), bool):
        state["dem_domain"] = value["dem_domain"]

    outlets = value.get("outlets")
    if isinstance(outlets, Mapping):
        state["outlets"] = {
            str(outlet_id): dict(record)
            for outlet_id, record in outlets.items()
            if isinstance(record, Mapping)
        }
    return state


def load_state(project_folder: PathInput) -> Dict[str, Any]:
    """Load state, returning defaults for missing or malformed content."""
    path = state_path(project_folder)
    try:
        with path.open("r", encoding="utf-8") as state_file:
            return _normalized_state(json.load(state_file))
    except (OSError, ValueError, TypeError):
        return default_state()


def _absolute_path(base: Path, value: PathInput) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def serialize_input_path(project_folder: PathInput, value: PathInput) -> str:
    """Serialize an input path relative to the outer project when possible."""
    if not str(value):
        return ""
    project = Path(project_folder).resolve()
    path = _absolute_path(project, value)
    try:
        return path.relative_to(project).as_posix()
    except ValueError:
        return str(path)


def resolve_input_path(project_folder: PathInput, value: PathInput) -> Path:
    """Resolve a stored input path relative to the outer project."""
    return _absolute_path(Path(project_folder).resolve(), value)


def serialize_output_path(project_folder: PathInput, value: PathInput) -> str:
    """Serialize a generated output relative to the plugin workspace."""
    if not str(value):
        return ""
    workspace = Path(workspace_folder(project_folder)).resolve()
    path = _absolute_path(workspace, value)
    try:
        return path.relative_to(workspace).as_posix()
    except ValueError as error:
        message = "Domain outputs must be inside the plugin workspace."
        raise ValueError(message) from error


def resolve_output_path(project_folder: PathInput, value: PathInput) -> Path:
    """Resolve a stored output path, rejecting paths outside the workspace."""
    workspace = Path(workspace_folder(project_folder)).resolve()
    path = _absolute_path(workspace, value)
    try:
        path.relative_to(workspace)
    except ValueError as error:
        message = "Domain outputs must be inside the plugin workspace."
        raise ValueError(message) from error
    return path


def _serialized_state(
    project_folder: PathInput,
    value: Mapping[str, Any],
) -> Dict[str, Any]:
    state = _normalized_state(value)
    for record in state["outlets"].values():
        if record.get("discharge_file"):
            record["discharge_file"] = serialize_input_path(
                project_folder, record["discharge_file"]
            )
        for key in ("mask_path", "vector_path"):
            if record.get(key):
                record[key] = serialize_output_path(project_folder, record[key])
    return state


def save_state(
    project_folder: PathInput,
    value: Mapping[str, Any],
) -> Path:
    """Atomically save normalized state and return its path."""
    path = state_path(project_folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    state = _serialized_state(project_folder, value)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8") as state_file:
            json.dump(state, state_file, indent=2, sort_keys=True)
            state_file.write("\n")
            state_file.flush()
            os.fsync(state_file.fileno())
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)
    return path


def _enabled(record: Mapping[str, Any], key: str, legacy_key: str) -> bool:
    return bool(record.get(key, record.get(legacy_key, False)))


def gauged_outlet_ids(state: Mapping[str, Any]) -> List[str]:
    """Return IDs of outlets marked as gauged, in stored order."""
    outlets = state.get("outlets", {})
    if not isinstance(outlets, Mapping):
        return []
    return [
        str(outlet_id)
        for outlet_id, record in outlets.items()
        if isinstance(record, Mapping)
        and _enabled(record, "is_gauged", "gauged")
    ]


def active_domain_records(state: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """Return active outlet domains plus the optional synthetic DEM domain."""
    records: List[Dict[str, Any]] = []
    outlets = state.get("outlets", {})
    if isinstance(outlets, Mapping):
        for outlet_id, record in outlets.items():
            if not isinstance(record, Mapping):
                continue
            if _enabled(record, "is_domain", "domain"):
                active = dict(record)
                active["outlet_id"] = str(outlet_id)
                active["is_dem_domain"] = False
                records.append(active)

    if state.get("dem_domain") is True:
        records.append(
            {
                "outlet_id": DEM_DOMAIN_ID,
                "is_domain": True,
                "is_dem_domain": True,
                "gauged_outlet_ids": gauged_outlet_ids(state),
            }
        )
    return records


def domain_count(state: Mapping[str, Any]) -> int:
    """Return the number of active outlet and DEM domains."""
    return len(active_domain_records(state))

"""Project-local inputs passed to the namelist GUI."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any

from ....project_layout import workspace_folder


def settings_path(project_folder: str | Path) -> Path:
    return Path(workspace_folder(project_folder)) / "nml-settings.json"


def load_settings(project_folder: str | Path) -> dict[str, Any]:
    path = settings_path(project_folder)
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        value = {}
    if not isinstance(value, dict):
        value = {}
    value.setdefault("version", 1)
    return value


def update_section(
    project_folder: str | Path,
    section: str,
    value: dict[str, Any] | list[Any],
) -> Path:
    state = load_settings(project_folder)
    state[section] = value
    return save_settings(project_folder, state)


def save_settings(
    project_folder: str | Path,
    state: dict[str, Any],
) -> Path:
    path = settings_path(project_folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(state)
    payload["version"] = 1
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return path


def relative_workspace_path(project_folder: str | Path, value: str | Path) -> str:
    return os.path.relpath(value, workspace_folder(project_folder)).replace("\\", "/")


def sync_domain_settings(project_folder: str | Path) -> Path:
    """Project saved domain and gauge records into the namelist handoff."""
    from .domain_state import (
        active_domain_records,
        gauge_records,
        load_state as load_domain_state,
    )

    domain_state = load_domain_state(project_folder)
    state = load_settings(project_folder)
    state["domain_definition"] = {
        "mode": str(domain_state.get("definition_mode", "") or ""),
        "dem_domain": bool(domain_state.get("dem_domain", False)),
    }
    gauges = gauge_records(domain_state)
    domains = []
    for record in active_domain_records(domain_state):
        if record.get("domain_id") is None:
            continue
        domain_id = int(record["domain_id"])
        outlet_id = str(record.get("outlet_id", ""))
        directory = str(record.get("domain_directory", "") or "")
        dem_path = str(record.get("dem_path", "") or "")
        domains.append(
            {
                "domain_id": domain_id,
                "name": "dem_extent" if record.get("is_dem_domain") else outlet_id,
                "outlet_id": outlet_id,
                "is_dem_domain": bool(record.get("is_dem_domain", False)),
                "directory": directory,
                "dem_path": dem_path,
                "data_path": "data/master/",
                "gauge_ids": [
                    int(gauge["gauge_id"])
                    for gauge in gauges
                    if gauge.get("gauge_id") is not None
                    and domain_id in gauge.get("domain_ids", [])
                ],
            }
        )
    state["domains"] = domains
    state["gauges"] = gauges
    return save_settings(project_folder, state)


__all__ = [
    "load_settings",
    "relative_workspace_path",
    "save_settings",
    "settings_path",
    "sync_domain_settings",
    "update_section",
]

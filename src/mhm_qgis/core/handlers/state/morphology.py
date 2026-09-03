# -*- coding: utf-8 -*-
"""The morphology half of `mhm_qgis_processing_state.json`.

Three writers share that file: `cache` records stage fingerprints,
`meteo_outputs` records prepared forcing, and this records prepared morphology
outputs, workflow status, the domain plan and the L0/L2 grid contract. Each
overlays only its own sections -- dumping a whole in-memory copy erases what the
others added since it was loaded, which silently disabled stage reuse once.

Functions take the session, so nothing here needs a dialog or a processor.
"""
from __future__ import annotations

import os

from ...utils.time import utc_timestamp
from ..file import json as jsonio
from ..store.paths import workspace_folder
from ..store.registry import available, key_for, register

STATE_FILENAME = "mhm_qgis_processing_state.json"

#: The only sections this writer owns in the shared state file.
OWNED_SECTIONS = ("version", "outputs", "workflows", "grid", "domains")

#: What an absent or unreadable state file looks like.
EMPTY = {"version": 1, "outputs": {}, "workflows": {}, "grid": {}}


def state_path(session):
    """Return the project-local processing state file, or None without a project."""
    if not session.project_folder:
        return None
    return os.path.join(workspace_folder(session.project_folder), STATE_FILENAME)


def load(session) -> dict:
    """Load the morphology registry onto the session and return it."""
    path = state_path(session)
    if not path or not os.path.exists(path):
        session.processing_state = dict(EMPTY)
        return session.processing_state

    state = jsonio.read(path)
    if not isinstance(state, dict):
        session.say("WARNING: Could not read processing state: not a JSON object.")
        session.processing_state = dict(EMPTY)
        return session.processing_state

    for key, default in EMPTY.items():
        state.setdefault(key, default if not isinstance(default, dict) else {})
    session.processing_state = state
    return state


def save(session) -> None:
    """Write the sections this registry owns, preserving the others.

    Atomic, via `jsonio.merge_sections`: this used to write the file in place,
    so a crash mid-write could corrupt state the other two writers share.
    """
    path = state_path(session)
    if not path:
        return
    state = session.processing_state or {}
    owned = {name: state[name] for name in OWNED_SECTIONS if name in state}
    try:
        jsonio.merge_sections(path, owned)
    except Exception as error:
        session.say(f"WARNING: Could not save processing state: {error}")


def output_key(session, path) -> str:
    """Return a stable registry key for an output path."""
    return key_for(session.project_folder, path)


def mark_prepared(session, path, name=None, loaded=False, algorithm=None) -> None:
    """Record that an output file has been prepared."""
    entry = register(session.project_folder, path,
                     name=name, loaded=loaded, algorithm=algorithm)
    if entry is not None:
        session.processing_state.setdefault("outputs", {})[entry["path"]] = entry


def is_prepared(session, path) -> bool:
    """Return True when an output is recorded and present on disk."""
    return available(session.project_folder, path)


def record_outputs(session, algorithm, params, result) -> None:
    """Record the file outputs a processing call declared."""
    if not result:
        return

    produced = [value for key, value in params.items()
                if key.upper().startswith("OUTPUT") and isinstance(value, str)]
    if isinstance(result, dict):
        produced += [value for key, value in result.items()
                     if key.upper().startswith("OUTPUT") and isinstance(value, str)]

    for path in produced:
        if path and path != "TEMPORARY_OUTPUT":
            mark_prepared(session, path, name=os.path.basename(path),
                          loaded=False, algorithm=algorithm)


def mark_workflow(session, workflow, status, message="", **metadata) -> None:
    """Record a project-local workflow status such as execute-all completion."""
    if not workflow:
        return

    timestamp = utc_timestamp()
    workflows = session.processing_state.setdefault("workflows", {})
    entry = workflows.get(workflow, {})
    entry.update({"status": status, "updated_at": timestamp})
    if status == "running":
        entry["started_at"] = timestamp
        entry.pop("completed_at", None)
        entry.pop("failed_at", None)
    elif status == "completed":
        entry["completed_at"] = timestamp
    elif status == "failed":
        entry["failed_at"] = timestamp

    if message:
        entry["message"] = message
    elif status == "running":
        entry.pop("message", None)

    for key, value in metadata.items():
        if value is not None:
            entry[key] = value

    workflows[workflow] = entry
    save(session)


def workflow_status(session, workflow) -> dict:
    """Return a saved workflow status entry."""
    return session.processing_state.get("workflows", {}).get(workflow, {})


def save_domain_plan(session, plan) -> list:
    """Record each domain's mask and target DEM for Morphology Setup."""
    session.processing_state["domains"] = [
        {key: entry[key] for key in
         ("domain_id", "outlet_id", "name", "mask", "directory", "dem_path")}
        for entry in plan
    ]
    save(session)
    return session.processing_state["domains"]


def saved_domain_plan(session) -> list:
    """Return the recorded domain plan."""
    plan = session.processing_state.get("domains")
    return list(plan) if isinstance(plan, list) else []


def save_grid_contract(session, l0_header, l2_header, multiplier) -> dict:
    """Record the validated L0/L2 headers and multiplier for later resumes."""
    from ....grid_resolution import validate_l0_l2_alignment

    ratio = int(multiplier)
    validate_l0_l2_alignment(l0_header, l2_header, ratio)
    session.processing_state["grid"] = {
        "l0_header": dict(l0_header),
        "l2_header": dict(l2_header),
        "l2_ratio_to_l0": ratio,
        "updated_at": utc_timestamp(),
    }
    save(session)
    return session.processing_state["grid"]


def saved_grid_contract(session):
    """Return the saved L0/L2 grid contract, or None when it is unusable."""
    from ....grid_resolution import validate_l0_l2_alignment

    grid = session.processing_state.get("grid") or {}
    l0_header, l2_header = grid.get("l0_header"), grid.get("l2_header")
    if not isinstance(l0_header, dict) or not isinstance(l2_header, dict):
        return None
    try:
        validate_l0_l2_alignment(l0_header, l2_header, int(grid.get("l2_ratio_to_l0")))
    except (TypeError, ValueError) as error:
        session.say(f"WARNING: Discarding the saved L0/L2 grid contract: {error}")
        return None
    return grid


__all__ = [
    "EMPTY", "OWNED_SECTIONS", "STATE_FILENAME",
    "is_prepared", "load", "mark_prepared", "mark_workflow", "output_key",
    "record_outputs", "save", "save_domain_plan", "save_grid_contract",
    "saved_domain_plan", "saved_grid_contract", "state_path", "workflow_status",
]

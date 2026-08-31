# -*- coding: utf-8 -*-
"""Where store finds the state file it shares with `core/handlers/state`.

Kept in its own module so `registry` never imports a state module directly --
that would close a loop, since the state modules need `store.paths` to locate
the workspace in the first place.
"""
from __future__ import annotations

from pathlib import Path

from .paths import workspace_folder

PROCESSING_STATE_FILENAME = "mhm_qgis_processing_state.json"


def processing_state_path(project_folder) -> Path:
    """Return the shared processing state path for one project."""
    return Path(workspace_folder(project_folder)) / PROCESSING_STATE_FILENAME

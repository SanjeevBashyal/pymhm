# -*- coding: utf-8 -*-
"""Persistent state for Configuration tab settings."""

from __future__ import annotations

import json
import os

from ..project_layout import (
    data_folder,
    geometry_folder,
    output_folder,
    restart_folder,
    z_temp_folder,
)
from ..time_utils import utc_timestamp


CONFIGURATION_STATE_FILENAME = "pymhm_configuration_state.json"


class ConfigurationStateStore:
    """Read and write all Configuration tab settings in one project JSON file."""

    def __init__(self, dialog):
        self.dialog = dialog

    def path(self):
        """Return the project-local configuration state file."""
        if not self.dialog.project_folder:
            return None
        return os.path.join(
            self.dialog.project_folder, CONFIGURATION_STATE_FILENAME)

    def empty(self):
        """Return an empty configuration state."""
        return {
            "version": 1,
            "mhm_version": "",
            "settings": {
                "mhm": {},
                "parameters": {},
                "outputs": {},
            },
        }

    def load(self):
        """Load configuration settings from disk."""
        path = self.path()
        if not path or not os.path.exists(path):
            return self.empty()

        try:
            with open(path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if not isinstance(state, dict):
                raise ValueError("Configuration state is not a JSON object.")
            state.setdefault("version", 1)
            state.setdefault("mhm_version", "")
            settings = state.setdefault("settings", {})
            settings.setdefault("mhm", {})
            settings.setdefault("parameters", {})
            settings.setdefault("outputs", {})
            return state
        except Exception as exc:
            self.dialog.log_message(
                f"WARNING: Could not read configuration state: {exc}")
            return self.empty()

    def save(self, state):
        """Save configuration settings to disk."""
        path = self.path()
        if not path:
            return

        state = dict(state or self.empty())
        state["version"] = 1
        state["updated_at"] = utc_timestamp()
        state["project_layout"] = self.project_layout()
        try:
            with open(path, "w", encoding="utf-8") as state_file:
                json.dump(state, state_file, indent=2, sort_keys=True)
        except Exception as exc:
            self.dialog.log_message(
                f"WARNING: Could not save configuration state: {exc}")

    def project_layout(self):
        """Return a serializable snapshot of the current project layout."""
        project_folder = self.dialog.project_folder
        if not project_folder:
            return {}
        return {
            "data_folder": data_folder(project_folder),
            "z_temp_folder": z_temp_folder(project_folder),
            "geometry_folder": geometry_folder(project_folder),
            "output_folder": output_folder(project_folder),
            "restart_folder": restart_folder(project_folder),
        }

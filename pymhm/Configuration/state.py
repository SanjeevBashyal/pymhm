# -*- coding: utf-8 -*-
"""Persistent state for Configuration tab settings."""

from __future__ import annotations

import json
import os
import tempfile
from typing import Any

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

    def __init__(self, dialog: Any) -> None:
        self.dialog = dialog

    def path(self) -> str | None:
        """Return the project-local configuration state file."""
        if not self.dialog.project_folder:
            return None
        return os.path.join(
            self.dialog.project_folder, CONFIGURATION_STATE_FILENAME)

    def empty(self) -> dict[str, Any]:
        """Return an empty configuration state."""
        return {
            "version": 1,
            "mhm_version": "",
        }

    def load(self) -> dict[str, Any]:
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
            settings = state.pop("settings", None)
            if isinstance(settings, dict):
                state["_legacy_namelist_settings"] = settings
            return state
        except Exception as exc:
            self.dialog.log_message(
                f"WARNING: Could not read configuration state: {exc}")
            return self.empty()

    def save(self, state: dict[str, Any]) -> None:
        """Save configuration settings to disk."""
        path = self.path()
        if not path:
            return

        state = {
            key: value
            for key, value in dict(state or self.empty()).items()
            if key not in ("settings", "_legacy_namelist_settings")
        }
        state["version"] = 1
        state["updated_at"] = utc_timestamp()
        state["project_layout"] = self.project_layout()
        try:
            directory = os.path.dirname(path) or "."
            os.makedirs(directory, exist_ok=True)
            descriptor, temporary = tempfile.mkstemp(
                prefix=".pymhm-state-", dir=directory)
            try:
                with os.fdopen(descriptor, "w", encoding="utf-8") as state_file:
                    json.dump(state, state_file, indent=2, sort_keys=True)
                    state_file.write("\n")
                os.replace(temporary, path)
            except Exception:
                try:
                    os.unlink(temporary)
                except OSError:
                    pass
                raise
        except Exception as exc:
            self.dialog.log_message(
                f"WARNING: Could not save configuration state: {exc}")

    def project_layout(self) -> dict[str, str]:
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

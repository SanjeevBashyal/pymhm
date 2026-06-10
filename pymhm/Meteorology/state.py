# -*- coding: utf-8 -*-
"""Project-local registry updates for meteorology outputs."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Union

from ..time_utils import utc_timestamp


STATE_FILENAME = "pymhm_processing_state.json"
PathInput = Union[str, Path]


class MeteorologyOutputState:
    """Small helper for recording prepared meteorology outputs."""

    def __init__(self, dialog: object) -> None:
        self.dialog = dialog

    def state_path(self) -> Path | None:
        if not self.dialog.project_folder:
            return None
        return Path(self.dialog.project_folder) / STATE_FILENAME

    def load(self) -> dict[str, object]:
        path = self.state_path()
        if not path or not path.exists():
            return {"version": 1, "outputs": {}}

        try:
            with path.open("r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if not isinstance(state, dict):
                raise ValueError("Processing state is not a JSON object.")
            state.setdefault("version", 1)
            state.setdefault("outputs", {})
            return state
        except Exception as e:
            self.dialog.log_message(
                f"WARNING: Could not read processing state: {e}")
            return {"version": 1, "outputs": {}}

    def save(self, state: dict[str, object]) -> None:
        path = self.state_path()
        if not path:
            return

        try:
            with path.open("w", encoding="utf-8") as state_file:
                json.dump(state, state_file, indent=2, sort_keys=True)
        except Exception as e:
            self.dialog.log_message(
                f"WARNING: Could not save processing state: {e}")

    def output_key(self, path: Path) -> str:
        try:
            return os.path.relpath(
                str(path), self.dialog.project_folder).replace("\\", "/")
        except ValueError:
            return str(path.resolve()).replace("\\", "/")

    def mark_prepared(self, path: PathInput, name: str | None = None) -> None:
        path = Path(path)
        if not path.exists():
            return

        state = self.load()
        key = self.output_key(path)
        entry = state.setdefault("outputs", {}).get(key, {})
        entry.update({
            "path": key,
            "absolute_path": str(path.resolve()),
            "exists": True,
            "loaded": False,
            "category": "meteorology",
            "updated_at": utc_timestamp(),
        })
        if name:
            entry["name"] = name

        state["outputs"][key] = entry
        self.save(state)

# -*- coding: utf-8 -*-
"""Project-local registry updates for meteorology outputs."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Union

from ..file import json as jsonio
from ....project_layout import workspace_folder
from ....time_utils import utc_timestamp


STATE_FILENAME = "mhm_qgis_processing_state.json"
PathInput = Union[str, Path]


class MeteorologyOutputState:
    """Small helper for recording prepared meteorology outputs."""

    def __init__(self, dialog: object) -> None:
        self.dialog = dialog

    def state_path(self) -> Path | None:
        if not self.dialog.project_folder:
            return None
        return Path(workspace_folder(self.dialog.project_folder)) / STATE_FILENAME

    #: The only sections this writer owns in the shared processing state.
    OWNED_SECTIONS = ("version", "outputs")

    def load(self) -> dict[str, object]:
        path = self.state_path()
        if not path:
            return {"version": 1, "outputs": {}}
        return jsonio.read_mapping(path, version=1, outputs={})

    def save(self, state: dict[str, object]) -> None:
        """Write back only this writer's sections of the shared state file.

        `cache` and the morphology processing state write their own sections
        into the same file. Dumping the whole in-memory copy erased whatever
        they had added since it was loaded, which silently disabled stage
        reuse -- so only the owned sections are overlaid, atomically.
        """
        path = self.state_path()
        if not path:
            return
        owned = {k: state[k] for k in self.OWNED_SECTIONS if k in state}
        try:
            jsonio.merge_sections(path, owned)
        except OSError as error:
            self.dialog.log_message(
                f"WARNING: Could not save processing state: {error}")

    def output_key(self, path: Path) -> str:
        try:
            return os.path.relpath(
                str(path),
                workspace_folder(self.dialog.project_folder),
            ).replace("\\", "/")
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

"""Tests for project-local state files in the plugin workspace."""
from __future__ import annotations

from pathlib import Path

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from pymhm.Meteorology.state import MeteorologyOutputState  # noqa: E402
from pymhm.Morphology.core.processing_state import ProcessingStateMixin  # noqa: E402
from pymhm.project_layout import workspace_folder  # noqa: E402
# isort: on


class _Dialog:
    def __init__(self, project_folder: Path) -> None:
        self.project_folder = str(project_folder)
        self.messages = []

    def log_message(self, message: str) -> None:
        self.messages.append(message)


class _MorphologyState(ProcessingStateMixin):
    def __init__(self, dialog: _Dialog) -> None:
        self.dialog = dialog
        self.processing_state_filename = "pymhm_processing_state.json"
        self.processing_state = {"version": 1, "outputs": {}, "workflows": {}}
        self.log_message = dialog.log_message


def test_morphology_state_is_saved_in_workspace(tmp_path: Path) -> None:
    dialog = _Dialog(tmp_path)
    state = _MorphologyState(dialog)
    expected = Path(workspace_folder(tmp_path)) / "pymhm_processing_state.json"

    assert Path(state.processing_state_path()) == expected
    state.save_processing_state()

    assert expected.is_file()


def test_meteorology_state_and_output_keys_use_workspace(tmp_path: Path) -> None:
    dialog = _Dialog(tmp_path)
    state = MeteorologyOutputState(dialog)
    workspace = Path(workspace_folder(tmp_path))
    output = workspace / "data" / "meteo" / "pre" / "pre.nc"

    assert state.state_path() == workspace / "pymhm_processing_state.json"
    assert state.output_key(output) == "data/meteo/pre/pre.nc"
    state.save({"version": 1, "outputs": {}})

    assert state.state_path().is_file()

"""Focused tests for the nml-tools configuration bridge."""
from __future__ import annotations

import os
import sys
import types
from pathlib import Path

import nml_tools.gui

from pymhm import configuration_processor as controller_module
from pymhm.configuration_processor import ConfigurationProcessor
from pymhm.dependency_bootstrap import (
    configure_qtpy_api,
    dependency_from_requirement_line,
)


class ComboBox:
    def __init__(self, text: str) -> None:
        self.text = text

    def currentText(self) -> str:
        return self.text


class Terminal:
    def __init__(self) -> None:
        self.calls = []

    def run_command(self, command, **kwargs):
        self.calls.append((command, kwargs))
        return True


class Dialog:
    def __init__(self, project_folder: Path, version: str = "6.0") -> None:
        self.project_folder = str(project_folder)
        self.comboBox_mHMversion = ComboBox(version)
        self.messages = []
        self.saved = 0
        self.terminal = Terminal()

    def log_message(self, message) -> None:
        self.messages.append(message)

    def save_input_state(self) -> None:
        self.saved += 1

    def open_project_terminal(self):
        return self.terminal


def test_edit_namelists_forwards_project_values_and_dimensions(
        tmp_path: Path, monkeypatch) -> None:
    dialog = Dialog(tmp_path, "5.13")
    calls = []
    prepared = []
    adapter = types.ModuleType("pymhm.vSpecific")
    adapter.build_initial_values = lambda version, owner: {
        "main": {"config_project": {"project_details": version}}
    }
    adapter.build_dimensions = lambda version, owner: {"n_domains": 2}
    monkeypatch.setitem(sys.modules, "pymhm.vSpecific", adapter)
    monkeypatch.setattr(
        controller_module,
        "ensure_project_structure",
        lambda path, version: prepared.append((path, version)),
    )
    monkeypatch.setattr(
        nml_tools.gui,
        "launch_gui",
        lambda **kwargs: calls.append(kwargs) or 0,
    )

    assert ConfigurationProcessor(dialog).edit_namelists() is True
    assert prepared == [(str(tmp_path), "5.13")]
    assert calls == [{
        "schemas_dir": Path(controller_module.plugin_root()) / "nml-schemas" / "v5.13",
        "output_dir": str(tmp_path),
        "initial_values": {
            "main": {"config_project": {"project_details": "5.13"}}
        },
        "initial_dimensions": {"n_domains": 2},
    }]


def test_run_mhm_uses_project_terminal(tmp_path: Path, monkeypatch) -> None:
    dialog = Dialog(tmp_path)
    monkeypatch.setattr(
        controller_module,
        "ensure_project_structure",
        lambda path, version: [],
    )

    assert ConfigurationProcessor(dialog).run_mhm() is True
    assert dialog.terminal.calls == [(
        "mhm",
        {"cwd": str(tmp_path), "show": True},
    )]


def test_version_change_prepares_project_and_saves_state(
        tmp_path: Path, monkeypatch) -> None:
    dialog = Dialog(tmp_path, "5.13")
    prepared = []
    monkeypatch.setattr(
        controller_module,
        "ensure_project_structure",
        lambda path, version: prepared.append((path, version)),
    )

    ConfigurationProcessor(dialog).handle_version_changed()

    assert prepared == [(str(tmp_path), "5.13")]
    assert dialog.saved == 1


def test_qtpy_is_forced_to_pyqt5(monkeypatch) -> None:
    monkeypatch.setenv("QT_API", "pyside6")

    assert configure_qtpy_api() == "pyqt5"
    assert os.environ["QT_API"] == "pyqt5"


def test_gui_dependency_import_names() -> None:
    nml_tools = dependency_from_requirement_line("nml-tools[gui]")
    assert nml_tools is not None
    assert nml_tools.import_name == "nml_tools"
    assert dependency_from_requirement_line("qtpy").import_name == "qtpy"
    assert dependency_from_requirement_line("guidata").import_name == "guidata"

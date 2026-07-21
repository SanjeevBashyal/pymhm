# -*- coding: utf-8 -*-
"""Open nml-tools for the selected mHM project and run mHM."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .dependency_bootstrap import configure_qtpy_api
from .project_layout import ensure_project_structure, plugin_root, version_key


class ConfigurationProcessor:
    """Bridge the pymhm dialog to nml-tools and project execution."""

    def __init__(self, dialog: Any) -> None:
        self.dialog = dialog
        self.log_message = dialog.log_message

    def selected_version(self) -> str:
        """Return the selected mHM version, defaulting to v6."""
        combo_box = getattr(self.dialog, "comboBox_mHMversion", None)
        if combo_box is None:
            return "6.0"
        return combo_box.currentText().strip() or "6.0"

    def handle_version_changed(self) -> None:
        """Prepare the selected project layout and persist the selection."""
        project_folder = self.dialog.project_folder
        if project_folder:
            ensure_project_structure(project_folder, self.selected_version())
        self.dialog.save_input_state()

    def edit_namelists(self) -> bool:
        """Open nml-tools with pymhm paths, values, and dimensions."""
        if not self.ensure_project_folder():
            return False

        try:
            from .vSpecific import build_dimensions, build_initial_values

            version = self.selected_version()
            project_folder = self.dialog.project_folder
            schemas_dir = Path(plugin_root()) / "nml-schemas" / version_key(version)
            ensure_project_structure(project_folder, version)
            configure_qtpy_api()

            from nml_tools.gui import launch_gui

            self.log_message(f"Opening namelist editor for mHM {version}...")
            launch_gui(
                schemas_dir=schemas_dir,
                output_dir=project_folder,
                initial_values=build_initial_values(version, self.dialog),
                initial_dimensions=build_dimensions(version, self.dialog),
            )
            return True
        except Exception as exc:
            self.report_error("Edit Namelists", exc)
            return False

    def run_mhm(self) -> bool:
        """Run mHM in the selected project folder through its terminal."""
        if not self.ensure_project_folder():
            return False

        project_folder = self.dialog.project_folder
        ensure_project_structure(project_folder, self.selected_version())
        terminal = self.dialog.open_project_terminal()
        if terminal is None:
            return False
        self.log_message(f"Running mHM in project directory: {project_folder}")
        return terminal.run_command("mhm", cwd=project_folder, show=True)

    def ensure_project_folder(self) -> bool:
        """Require a selected project folder."""
        if self.dialog.project_folder:
            return True
        self._message_box().warning(
            self.dialog,
            "Project Folder Required",
            "Select a project folder before editing namelists.",
        )
        return False

    def report_error(self, title: str, exc: Exception) -> None:
        """Log and show a configuration error."""
        self.log_message(f"ERROR: {title} failed. Details: {exc}")
        self._message_box().critical(self.dialog, "Configuration Error", str(exc))

    @staticmethod
    def _message_box() -> Any:
        from qgis.PyQt.QtWidgets import QMessageBox

        return QMessageBox


__all__ = ["ConfigurationProcessor"]

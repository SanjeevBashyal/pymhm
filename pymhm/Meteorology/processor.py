# -*- coding: utf-8 -*-
"""Meteorology button actions for pymhm."""
from __future__ import annotations

from typing import Any

from qgis.PyQt.QtWidgets import QMessageBox

from .era5_tools import load_era5_mhm_tools
from .exceptions import MeteorologyInputError, MeteorologyToolImportError
from .inputs import collect_run_inputs
from .inventory import log_inventory
from .paths import expected_meteo_outputs
from .results import log_result_summary, record_forcing_outputs
from .state import MeteorologyOutputState
from ..grid_resolution import save_meteo_grid_metadata


class MeteorologyProcessor:
    """Handle meteorology forcing preparation from the pymhm dialog."""

    def __init__(self, dialog: Any) -> None:
        self.dialog = dialog
        self.log_message = dialog.log_message
        self.state = MeteorologyOutputState(dialog)

    def load_project_state(self) -> None:
        """Record and report existing final meteorology forcing files."""
        if not self.dialog.project_folder:
            return

        found = []
        for variable, paths in expected_meteo_outputs(
                self.dialog.project_folder).items():
            if paths["netcdf"].exists() and paths["header"].exists():
                self.state.mark_prepared(paths["netcdf"], paths["netcdf"].name)
                self.state.mark_prepared(paths["header"], paths["header"].name)
                found.append(variable)

        if found:
            self.log_message(
                "Existing meteorology forcing files found: "
                + ", ".join(found))
        if hasattr(self.dialog, "update_l2_resolution_from_metadata"):
            self.dialog.update_l2_resolution_from_metadata()

    def process_meteo_forcing(self) -> bool:
        """Prepare mHM meteo forcing NetCDF files from ERA5-Land inputs."""
        self.log_message("\n--- Preparing meteorology forcing files ---")

        try:
            run_inputs = collect_run_inputs(self.dialog)
        except MeteorologyInputError as e:
            self._show_input_error(e)
            return False

        self.log_message(f"ERA5-Land input folder: {run_inputs.nc_folder}")
        self.log_message(f"mHM meteorology output folder: {run_inputs.output_root}")

        try:
            tools = load_era5_mhm_tools()
        except MeteorologyToolImportError as e:
            QMessageBox.critical(
                self.dialog,
                "Meteorology Error",
                str(e))
            self.log_message(f"ERROR: Could not import ERA5-Land tools: {e}")
            return False

        missing_dependency_error = tools.missing_dependency_error
        inventory = tools.inspect_folder(run_inputs.nc_folder)
        log_inventory(inventory, self.log_message)

        try:
            l2_grid = self.dialog.prepare_meteo_l2_grid(run_inputs.nc_folder)
            result = tools.process_to_mhm(
                nc_folder=run_inputs.nc_folder,
                output_root=run_inputs.output_root,
                bounds=l2_grid.get("bounds", run_inputs.crop_bounds),
                target_lat=l2_grid["lat"],
                target_lon=l2_grid["lon"],
                target_header=l2_grid["header"],
                skip_existing=False,
                log=self.log_message,
            )
        except missing_dependency_error as e:
            self.log_message(f"ERROR: {e}")
            QMessageBox.warning(
                self.dialog,
                "Missing Dependency",
                f"{e}\n\nInstall the missing packages in the QGIS Python environment.")
            return False
        except Exception as e:
            import traceback
            details = traceback.format_exc()
            self.log_message(f"ERROR: Meteorology processing failed: {e}\n{details}")
            QMessageBox.critical(
                self.dialog,
                "Meteorology Error",
                f"Meteorology processing failed.\n{e}")
            return False

        record_forcing_outputs(self.state, result)
        log_result_summary(result, self.log_message)
        save_meteo_grid_metadata(
            run_inputs.project_folder,
            l2_grid["metadata"],
        )
        if hasattr(self.dialog, "set_meteo_l2_grid_metadata"):
            self.dialog.set_meteo_l2_grid_metadata(l2_grid["metadata"])

        QMessageBox.information(
            self.dialog,
            "Success",
            "Meteorology forcing files prepared successfully.")
        return True

    def _show_input_error(self, error: MeteorologyInputError) -> None:
        """Show validated input errors with the appropriate dialog severity."""
        if error.severity == "critical":
            QMessageBox.critical(self.dialog, error.title, error.message)
        else:
            QMessageBox.warning(self.dialog, error.title, error.message)

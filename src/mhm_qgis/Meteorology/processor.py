# -*- coding: utf-8 -*-
"""Meteorology preparation orchestration for the mhm_qgis dialog."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import traceback
from typing import Any

from qgis.PyQt.QtWidgets import QMessageBox

from .forcing import (MeteoFolderSpec, TargetGrid,
                      process_meteo_inputs)
from .paths import expected_meteo_outputs, meteo_mask_path, meteo_output_root
from .results import log_result_summary, record_forcing_outputs
from .reuse import required_meteo_variables, stale_meteo_variables
from .state import MeteorologyOutputState
from .mask import write_meteo_mask
from ..grid_resolution import save_meteo_grid_metadata
from ..project_layout import is_v6


@dataclass(frozen=True)
class MeteorologyRun:
    """Immutable inputs prepared on the UI thread for one meteo run."""

    project_folder: Path
    precipitation: MeteoFolderSpec
    temperature: MeteoFolderSpec
    pet: MeteoFolderSpec | None
    target_grid: TargetGrid
    grid_metadata: dict


class MeteorologyProcessor:
    """Compile, resample, crop, and write meteorology forcing."""

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
            header = paths["header"]
            if not paths["netcdf"].exists():
                continue
            if header is not None and not header.exists():
                continue
            self.state.mark_prepared(paths["netcdf"], paths["netcdf"].name)
            if header is not None:
                self.state.mark_prepared(header, header.name)
            found.append(variable)

        if found:
            self.log_message(
                "Existing meteorology forcing files found: "
                + ", ".join(found))
        if hasattr(self.dialog, "update_l2_resolution_from_metadata"):
            self.dialog.update_l2_resolution_from_metadata()

    def reuse_prepared_meteo_forcing(self, run: MeteorologyRun) -> bool:
        """Return True when recorded forcing already matches the target L2 grid."""
        variables = required_meteo_variables(run.pet is not None)
        stale = stale_meteo_variables(
            run.project_folder,
            self.state.load().get("outputs", {}),
            run.target_grid.header,
            variables,
        )
        if stale:
            for variable, reason in sorted(stale.items()):
                self.log_message(f"{variable}: needs preparing ({reason}).")
            return False

        self.log_message(
            "Meteorology forcing is already prepared on this L2 grid for "
            f"{', '.join(variables)}. Skipping meteorology preparation."
        )
        self.log_message(
            "Change the L2 grid or delete the files under "
            f"{meteo_output_root(run.project_folder)} to force a rebuild."
        )
        save_meteo_grid_metadata(run.project_folder, run.grid_metadata)
        return True

    def process_meteo_forcing(
            self,
            run: MeteorologyRun,
            show_dialog: bool = True) -> bool:
        """Prepare all selected meteorology forcing files."""
        self.log_message("\n--- Preparing meteorology forcing files ---")
        self.log_message(
            f"Meteorology output folder: "
            f"{meteo_output_root(run.project_folder)}"
        )
        if self.reuse_prepared_meteo_forcing(run):
            return True

        try:
            flat = is_v6(run.project_folder)
            result = process_meteo_inputs(
                precipitation=run.precipitation,
                temperature=run.temperature,
                pet=run.pet,
                output_root=meteo_output_root(run.project_folder),
                target_grid=run.target_grid,
                flat_layout=flat,
                log=self.log_message,
            )
            if flat:
                # v6 names a meteorology mask in config_input.
                mask = write_meteo_mask(
                    meteo_mask_path(run.project_folder), run.target_grid)
                self.state.mark_prepared(mask, mask.name)
                self.log_message(f"Meteorology mask written: {mask}")
        except Exception as error:
            self.log_message(
                f"ERROR: Meteorology processing failed: {error}\n"
                f"{traceback.format_exc()}"
            )
            if show_dialog:
                QMessageBox.critical(
                    self.dialog,
                    "Meteorology Error",
                    f"Meteorology processing failed.\n{error}",
                )
            return False

        record_forcing_outputs(self.state, result)
        log_result_summary(result, self.log_message)
        save_meteo_grid_metadata(run.project_folder, run.grid_metadata)
        if show_dialog:
            self.dialog.set_meteo_l2_grid_metadata(run.grid_metadata)
            QMessageBox.information(
                self.dialog,
                "Success",
                "Meteorology forcing files prepared successfully.",
            )
        return True

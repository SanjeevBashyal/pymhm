# -*- coding: utf-8 -*-
"""Hydrology page discharge-table assignment workflow."""
from __future__ import annotations

from ..common import QMessageBox
from ..core.base import BaseProcessingMixin
from .discharge_dialog import DischargeTableAssignmentDialog
from .discharge_writer import write_streamflow_observation
from .observation_paths import streamflow_observation_folder
from .outlets import StationIdError, station_ids_from_layer


class DischargeAssignmentMixin(BaseProcessingMixin):
    """Assign selected discharge table layers to gauged outlets."""

    def assign_discharge_tables(self) -> None:
        """Open the discharge table assignment dialog and write mHM files."""
        self.log_message("\n--- Assigning discharge tables to gauged outlets ---")

        if not self.dialog.project_folder:
            QMessageBox.warning(
                self.dialog,
                "Missing Input",
                "Please select a project folder before assigning discharge tables.")
            return

        pour_points_layer = self.dialog.mMapLayerComboBox_pour_points.currentLayer()
        try:
            station_ids = station_ids_from_layer(pour_points_layer)
        except StationIdError as e:
            QMessageBox.critical(self.dialog, "Missing STATION_ID", str(e))
            self.log_message(f"ERROR: {e}")
            return

        if not station_ids:
            QMessageBox.warning(
                self.dialog,
                "No Gauged Outlets",
                "The selected pour points layer does not contain any features.")
            return

        dialog = DischargeTableAssignmentDialog(station_ids, self.dialog)
        if dialog.exec_() != dialog.Accepted:
            self.log_message("Discharge table assignment cancelled.")
            return

        selected_layers = dialog.selected_layers()
        output_folder = streamflow_observation_folder(self.dialog.project_folder)
        self.log_message(f"Streamflow observation output folder: {output_folder}")

        written = []
        for station_id in station_ids:
            layer = selected_layers.get(station_id)
            if not layer or not layer.isValid():
                QMessageBox.warning(
                    self.dialog,
                    "Missing Discharge Table",
                    f"Please select a valid discharge table for STATION_ID {station_id}.")
                return

            self.log_message(
                f"STATION_ID {station_id}: discharge table {layer.name()} | {layer.source()}")
            try:
                output_file = write_streamflow_observation(
                    layer,
                    station_id,
                    output_folder,
                )
            except Exception as e:
                QMessageBox.critical(
                    self.dialog,
                    "Discharge Table Error",
                    f"Could not write streamflow observations for STATION_ID {station_id}.\n{e}")
                self.log_message(
                    f"ERROR: Could not write streamflow observations for {station_id}: {e}")
                return

            self.mark_output_prepared(
                str(output_file),
                name=output_file.name,
                loaded=False,
            )
            written.append(output_file)
            self.log_message(f"Written streamflow observations: {output_file}")

        self.log_message(
            f"Discharge table assignment completed for {len(written)} gauged outlet(s).")
        QMessageBox.information(
            self.dialog,
            "Success",
            f"Prepared {len(written)} streamflow observation file(s).")

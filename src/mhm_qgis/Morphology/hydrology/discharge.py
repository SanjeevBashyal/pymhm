# -*- coding: utf-8 -*-
"""Hydrology page discharge-table assignment workflow."""
from __future__ import annotations

from ..common import QMessageBox
from ..core.base import BaseProcessingMixin
from .discharge_dialog import DischargeTableAssignmentDialog
from .outlets import (
    StationIdError,
    selected_outlet_id_field,
    station_ids_from_layer,
)
from ..watershed.domain_state import (
    active_domain_records,
    load_state,
    save_state,
    state_path,
)
from ..watershed.domain_workflow import DomainWorkflow


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

        pour_points_layer = self.dialog.input_combo("pour_points").currentLayer()
        try:
            all_station_ids = station_ids_from_layer(
                pour_points_layer,
                field_name=selected_outlet_id_field(self.dialog),
            )
            station_ids = all_station_ids
        except StationIdError as e:
            QMessageBox.critical(self.dialog, "Invalid Outlet IDs", str(e))
            self.log_message(f"ERROR: {e}")
            return

        if not station_ids:
            QMessageBox.warning(
                self.dialog,
                "No Pour Points",
                "The selected pour-point layer does not contain any features.")
            return

        existing_state = (
            load_state(self.dialog.project_folder)
            if state_path(self.dialog.project_folder).is_file()
            else None
        )
        dialog = DischargeTableAssignmentDialog(
            station_ids,
            self.dialog,
            initial_records=(existing_state or {}).get("outlets", {}),
        )
        if dialog.exec_() != dialog.Accepted:
            self.log_message("Discharge table assignment cancelled.")
            return

        assignments = dialog.selected_assignments()
        workflow = DomainWorkflow(
            self.dialog,
            self,
            pour_points_layer,
            selected_outlet_id_field(self.dialog),
            station_ids,
        )
        try:
            prepared = workflow.validate_gauge_assignments(assignments)
            if not prepared:
                raise ValueError("Select at least one gauged outlet.")
            state = None
            previous = set()
            if existing_state is not None:
                previous = set(
                    outlet_id
                    for outlet_id, record in existing_state["outlets"].items()
                    if record.get("is_gauged", record.get("gauged", False))
                )
                state = workflow.load_synced_state(
                    existing_state.get("definition_mode", ""),
                    existing_state.get("dem_domain", False),
                )
                workflow.apply_assignment_records(state, assignments, prepared)
                workflow.validate_unique_state_gauge_ids(state)
                if active_domain_records(state):
                    workflow.update_gauge_domain_ids(state)
                else:
                    for record in state["outlets"].values():
                        if record.get("is_gauged", False):
                            record["domain_ids"] = []

            workflow.write_gauges(prepared)
            if state is not None:
                save_state(self.dialog.project_folder, state)
                workflow.remove_deselected_gauges(previous, state)
        except (StationIdError, ValueError, RuntimeError) as error:
            QMessageBox.critical(
                self.dialog,
                "Discharge Table Error",
                str(error),
            )
            self.log_message(f"ERROR: {error}")
            return

        for gauge in prepared.values():
            self.log_message(
                f"Written streamflow observations: {gauge.output_path}"
            )

        self.log_message(
            "Discharge table assignment completed for "
            f"{len(prepared)} gauged outlet(s)."
        )
        QMessageBox.information(
            self.dialog,
            "Success",
            f"Prepared {len(prepared)} streamflow observation file(s).",
        )

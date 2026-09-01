# -*- coding: utf-8 -*-
"""Widget state for `domain_delineator_dialog.ui`.

The embedded map canvas, its map tool and the vertex marker stay on the dialog:
a `QgsMapToolEmitPoint` belongs to the widget that owns the canvas, and
`qt/bindings/domain_delineator.py` already draws that line.
"""
from __future__ import annotations

import os
from pathlib import Path

import copy

try:
    from qgis.core import QgsVectorLayer
    from qgis.PyQt.QtWidgets import QFileDialog, QMessageBox
except ImportError:
    from ...standalone import install

    install(force=True)
    from qgis.core import QgsVectorLayer
    from qgis.PyQt.QtWidgets import QFileDialog, QMessageBox

from ..dialogs.input_selection import scan_project_inputs
from ..dialogs.discharge_assignment import OutletAssignment
from ...Morphology.hydrology.outlets import StationIdError

def _load_outlet(dialog, outlet_id):
    dialog.current_outlet_id = str(outlet_id or "")
    if not dialog.current_outlet_id:
        dialog._outlet_marker.hide()
        return

    dialog._preview_result = None
    record = dialog._draft_state["outlets"].get(dialog.current_outlet_id, {})
    dialog._channel_layer = QgsVectorLayer(
        dialog.processor.channel_network_vector_path,
        "Channel network",
        "ogr",
    )
    dialog.label_outletIDValue.setText(dialog.current_outlet_id)
    dialog.checkBox_isGaugedOutlet.blockSignals(True)
    dialog.checkBox_isDomainOutlet.blockSignals(True)
    dialog.checkBox_isGaugedOutlet.setChecked(
        bool(record.get("is_gauged", record.get("gauged", False))))
    dialog.checkBox_isDomainOutlet.setChecked(
        bool(record.get("is_domain", record.get("domain", False))))
    dialog.checkBox_isGaugedOutlet.blockSignals(False)
    dialog.checkBox_isDomainOutlet.blockSignals(False)
    dialog._set_discharge_enabled(
        dialog.checkBox_isGaugedOutlet.isChecked())
    dialog.widget_domainControls.setEnabled(True)
    try:
        threshold = max(1, int(record.get("threshold_cells", 1) or 1))
    except (TypeError, ValueError):
        threshold = 1
    dialog.spinBox_channelThreshold.setValue(threshold)

    discharge = dialog._resolved_input(record.get("discharge_file"))
    dialog._populate_discharge_files(discharge)
    area = record.get("catchment_area_m2")
    dialog._show_area(area)
    dialog._show_picked_coordinates(record.get("picked"))
    dialog._show_saved_watershed(record)
    dialog._zoom_to_outlet()


def _set_discharge_enabled(dialog, enabled):
    dialog.comboBox_dischargeFile.setEnabled(bool(enabled))
    dialog.pushButton_browseDischargeFile.setEnabled(bool(enabled))


def _populate_discharge_files(dialog, current_path=""):
    current = dialog._normal_path(current_path) if current_path else ""
    excluded = dialog._used_discharge_paths(dialog.current_outlet_id)
    dialog.comboBox_dischargeFile.blockSignals(True)
    dialog.comboBox_dischargeFile.clear()
    dialog.comboBox_dischargeFile.addItem("", "")

    available = {}
    for item in scan_project_inputs(dialog.project_folder, "lookup"):
        path = dialog._normal_path(item.data["path"])
        if path in excluded:
            continue
        available[path] = item.label
    if current and current not in excluded and os.path.isfile(current):
        available.setdefault(current, dialog._path_label(current))

    selected_index = 0
    for path, label in sorted(
        available.items(), key=lambda item: item[1].casefold()
    ):
        dialog.comboBox_dischargeFile.addItem(label, path)
        if path == current:
            selected_index = dialog.comboBox_dischargeFile.count() - 1
    dialog.comboBox_dischargeFile.setCurrentIndex(selected_index)
    dialog.comboBox_dischargeFile.blockSignals(False)


def _path_label(dialog, path):
    try:
        return Path(path).resolve().relative_to(
            Path(dialog.project_folder).resolve()).as_posix()
    except ValueError:
        return str(Path(path).resolve())


def _browse_discharge(dialog):
    path, _ = QFileDialog.getOpenFileName(
        dialog,
        "Select discharge data",
        dialog.project_folder,
        "Discharge data (*.csv *.txt)",
    )
    if not path:
        return
    if Path(path).suffix.lower() not in {".csv", ".txt"}:
        QMessageBox.warning(
            dialog, "Invalid file", "Select a CSV or TXT discharge file.")
        return

    normalized = dialog._normal_path(path)
    if normalized in dialog._used_discharge_paths(dialog.current_outlet_id):
        QMessageBox.warning(
            dialog,
            "File already selected",
            "This file is already selected elsewhere in the plugin.",
        )
        return
    dialog._select_discharge_path(normalized)


def _select_discharge_path(dialog, path):
    for index in range(dialog.comboBox_dischargeFile.count()):
        if dialog.comboBox_dischargeFile.itemData(index) == path:
            dialog.comboBox_dischargeFile.setCurrentIndex(index)
            return
    dialog.comboBox_dischargeFile.addItem(dialog._path_label(path), path)
    dialog.comboBox_dischargeFile.setCurrentIndex(
        dialog.comboBox_dischargeFile.count() - 1)


def _current_assignment(dialog):
    outlet_id = dialog.current_outlet_id
    is_gauged = dialog.checkBox_isGaugedOutlet.isChecked()
    discharge_path = str(dialog.comboBox_dischargeFile.currentData() or "")
    if is_gauged and not discharge_path:
        raise ValueError(
            f"Select a discharge CSV or TXT file for gauge {outlet_id}."
        )
    return OutletAssignment(
        outlet_id=outlet_id,
        is_gauge=is_gauged,
        is_domain=dialog.checkBox_isDomainOutlet.isChecked(),
        discharge_layer=(
            QgsVectorLayer(discharge_path, Path(discharge_path).name, "ogr")
            if is_gauged and discharge_path
            else None
        ),
    )


def _next_outlet(dialog):
    try:
        dialog._stage_current_outlet()
    except (StationIdError, ValueError, RuntimeError) as error:
        dialog._show_error("Next Outlet", error)
        return
    row = dialog.listWidget_outlets.currentRow()
    if row + 1 < dialog.listWidget_outlets.count():
        dialog.listWidget_outlets.setCurrentRow(row + 1)


def _domain_state_saved(dialog, proposed_state):
    dialog.state = proposed_state
    dialog._draft_state = copy.deepcopy(proposed_state)
    dialog.main_dialog.save_input_state()
    dialog.processor.update_gauged_outlet_count()
    dialog._show_saved_watershed(
        proposed_state["outlets"].get(dialog.current_outlet_id, {})
    )
    QMessageBox.information(
        dialog, "Domain Delineator", "Saved all pour-point inputs."
    )


def _show_area(dialog, value):
    try:
        area = float(value)
    except (TypeError, ValueError):
        dialog.label_catchmentAreaValue.setText("-")
        return
    dialog.label_catchmentAreaValue.setText(
        f"{area / 1_000_000.0:.3f} km²")


def _show_picked_coordinates(dialog, picked):
    if not isinstance(picked, dict):
        dialog.label_pourPointUpdatedCoordinatesValue.setText("-")
        return
    try:
        x = float(picked["x"])
        y = float(picked["y"])
    except (KeyError, TypeError, ValueError):
        dialog.label_pourPointUpdatedCoordinatesValue.setText("-")
        return
    authid = str(picked.get("crs", "") or "")
    suffix = f" ({authid})" if authid else ""
    dialog.label_pourPointUpdatedCoordinatesValue.setText(
        f"{x:.3f}, {y:.3f}{suffix}"
    )


def _show_error(dialog, title, error):
    dialog.main_dialog.log_message(f"ERROR: {title}: {error}")
    QMessageBox.critical(dialog, title, str(error))

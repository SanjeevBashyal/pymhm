# -*- coding: utf-8 -*-
"""Per-outlet gauge and domain assignment dialogs."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from qgis.PyQt.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QLabel,
    QMessageBox,
    QPushButton,
)
from qgis.core import QgsMapLayer, QgsProject, QgsVectorLayer

from ...ui.pyui.ui_discharge_table_assignment_dialog import (
    Ui_DischargeTableAssignmentDialog,
)
from ...ui.pyui.ui_domain_and_discharge_table_assignment_dialog import (
    Ui_DischargeTableAssignmentDialog as Ui_DomainAndDischargeTableAssignmentDialog,
)


@dataclass(frozen=True)
class OutletAssignment:
    """One outlet row returned by an assignment dialog."""

    outlet_id: str
    is_gauge: bool
    is_domain: bool = False
    discharge_layer: object | None = None

    @property
    def discharge_source(self) -> str:
        """Return the selected layer source, if one is available."""
        source = getattr(self.discharge_layer, "source", None)
        return str(source() or "") if callable(source) else ""


@dataclass
class _AssignmentRow:
    outlet_id: str
    id_value: QLabel
    gauge: QCheckBox
    discharge: QComboBox
    browse: QPushButton
    domain: QCheckBox | None = None


class _OutletAssignmentRows:
    """Populate a Designer row template once per outlet."""

    include_domain = False
    default_is_gauge = False

    def _setup_assignment_rows(
        self,
        outlet_ids: list[str],
        initial_records: Mapping[str, Mapping] | None = None,
    ) -> None:
        self._rows: list[_AssignmentRow] = []
        self._initial_records = initial_records or {}
        self._remove_layout_spacers()

        for row_index, outlet_id in enumerate(outlet_ids):
            row = (
                self._designer_row(str(outlet_id))
                if row_index == 0
                else self._create_row(row_index, str(outlet_id))
            )
            self._configure_row(row)
            self._rows.append(row)

        if not outlet_ids:
            for widget in self._designer_row_widgets():
                widget.hide()

    def _remove_layout_spacers(self) -> None:
        for index in range(self.rows_layout.count() - 1, -1, -1):
            item = self.rows_layout.itemAt(index)
            if item is not None and item.spacerItem() is not None:
                self.rows_layout.takeAt(index)

    def _designer_row(self, outlet_id: str) -> _AssignmentRow:
        return _AssignmentRow(
            outlet_id=outlet_id,
            id_value=self.label_pourPointIDValue,
            domain=(self.checkBox_isDomain if self.include_domain else None),
            gauge=self.checkBox_isGauge,
            discharge=self.comboBox_dischargeTableInput,
            browse=self.pushButton_browseDischargeTable1,
        )

    def _designer_row_widgets(self) -> tuple[object, ...]:
        widgets = [
            self.label_pourPointID,
            self.label_pourPointIDValue,
            self.checkBox_isGauge,
            self.comboBox_dischargeTableInput,
            self.pushButton_browseDischargeTable1,
        ]
        if self.include_domain:
            widgets.append(self.checkBox_isDomain)
        return tuple(widgets)

    def _create_row(self, row_index: int, outlet_id: str) -> _AssignmentRow:
        suffix = row_index + 1
        label = QLabel("Pour Point ID:", self.rows_widget)
        label.setObjectName(f"label_pourPointID{suffix}")
        id_value = QLabel(self.rows_widget)
        id_value.setObjectName(f"label_pourPointIDValue{suffix}")
        gauge = QCheckBox("Is Gauge?", self.rows_widget)
        gauge.setObjectName(f"checkBox_isGauge{suffix}")
        discharge = QComboBox(self.rows_widget)
        discharge.setObjectName(f"comboBox_dischargeTableInput{suffix}")
        browse = QPushButton("...", self.rows_widget)
        browse.setObjectName(f"pushButton_browseDischargeTable{suffix}")
        browse.setMaximumWidth(40)

        column = 0
        self.rows_layout.addWidget(label, row_index, column)
        self.rows_layout.addWidget(id_value, row_index, column + 1)
        domain = None
        if self.include_domain:
            domain = QCheckBox("Is Domain?", self.rows_widget)
            domain.setObjectName(f"checkBox_isDomain{suffix}")
            self.rows_layout.addWidget(domain, row_index, column + 2)
            column += 1
        self.rows_layout.addWidget(gauge, row_index, column + 2)
        self.rows_layout.addWidget(discharge, row_index, column + 3)
        self.rows_layout.addWidget(browse, row_index, column + 4)
        return _AssignmentRow(
            outlet_id=outlet_id,
            id_value=id_value,
            domain=domain,
            gauge=gauge,
            discharge=discharge,
            browse=browse,
        )

    def _configure_row(self, row: _AssignmentRow) -> None:
        record = self._initial_records.get(row.outlet_id, {})
        row.id_value.setText(row.outlet_id)
        self._populate_discharge_layers(row.discharge)

        is_gauge = bool(
            record.get(
                "is_gauged",
                record.get("gauged", self.default_is_gauge),
            )
        )
        row.gauge.setChecked(is_gauge)
        if row.domain is not None:
            row.domain.setChecked(
                bool(record.get("is_domain", record.get("domain", False)))
            )

        source = str(record.get("discharge_file", "") or "")
        if source:
            self._select_discharge_source(row.discharge, source)
        row.gauge.toggled.connect(
            lambda enabled, current=row: self._set_discharge_enabled(
                current, enabled
            )
        )
        row.browse.clicked.connect(
            lambda checked=False, current=row: self._browse_discharge(current)
        )
        self._set_discharge_enabled(row, is_gauge)

    def _populate_discharge_layers(self, combo: QComboBox) -> None:
        combo.clear()
        combo.addItem("", None)
        try:
            layers = QgsProject.instance().mapLayers().values()
        except Exception:
            layers = ()
        for layer in layers:
            try:
                if layer.isValid() and layer.type() == QgsMapLayer.VectorLayer:
                    combo.addItem(layer.name(), layer)
            except Exception:
                continue

    def _select_discharge_source(self, combo: QComboBox, source: str) -> None:
        resolved = self._resolve_input_source(source)
        target = self._normal_source(resolved)
        for index in range(combo.count()):
            layer = combo.itemData(index)
            if layer is not None and self._normal_source(layer.source()) == target:
                combo.setCurrentIndex(index)
                return
        if os.path.isfile(resolved):
            layer = QgsVectorLayer(resolved, Path(resolved).name, "ogr")
            if layer.isValid():
                combo.addItem(layer.name(), layer)
                combo.setCurrentIndex(combo.count() - 1)

    def _resolve_input_source(self, source: str) -> str:
        path = str(source).split("|", 1)[0]
        project_folder = getattr(self.parent(), "project_folder", None)
        if path and project_folder and not os.path.isabs(path):
            path = os.path.join(str(project_folder), path)
        return os.path.abspath(path) if path else ""

    @staticmethod
    def _normal_source(source: str) -> str:
        path = str(source or "").split("|", 1)[0].split("?", 1)[0]
        return os.path.normcase(os.path.abspath(path)) if path else ""

    @staticmethod
    def _set_discharge_enabled(row: _AssignmentRow, enabled: bool) -> None:
        row.discharge.setEnabled(bool(enabled))
        row.browse.setEnabled(bool(enabled))

    def _browse_discharge(self, row: _AssignmentRow) -> None:
        project_folder = getattr(self.parent(), "project_folder", "") or ""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select discharge data",
            str(project_folder),
            "Discharge data (*.csv *.txt);;All files (*)",
        )
        if not path:
            return
        layer = QgsVectorLayer(path, Path(path).name, "ogr")
        if not layer.isValid():
            QMessageBox.warning(
                self,
                "Invalid discharge table",
                f"Could not read the selected table:\n{path}",
            )
            return
        target = self._normal_source(path)
        for index in range(row.discharge.count()):
            candidate = row.discharge.itemData(index)
            if (
                candidate is not None
                and self._normal_source(candidate.source()) == target
            ):
                row.discharge.setCurrentIndex(index)
                return
        row.discharge.addItem(layer.name(), layer)
        row.discharge.setCurrentIndex(row.discharge.count() - 1)

    def selected_assignments(self) -> list[OutletAssignment]:
        """Return one typed assignment for every displayed outlet."""
        return [
            OutletAssignment(
                outlet_id=row.outlet_id,
                is_gauge=row.gauge.isChecked(),
                is_domain=(row.domain.isChecked() if row.domain else False),
                discharge_layer=(
                    row.discharge.currentData() if row.gauge.isChecked() else None
                ),
            )
            for row in self._rows
        ]

    def selected_layers(self) -> dict[str, object | None]:
        """Backward-compatible selected-layer mapping keyed by outlet ID."""
        return {
            record.outlet_id: record.discharge_layer
            for record in self.selected_assignments()
        }


class DischargeTableAssignmentDialog(
    QDialog,
    Ui_DischargeTableAssignmentDialog,
    _OutletAssignmentRows,
):
    """Assign optional discharge tables to all supplied outlets."""

    default_is_gauge = True

    def __init__(
        self,
        station_ids: list[str],
        parent=None,
        initial_records: Mapping[str, Mapping] | None = None,
    ):
        super().__init__(parent)
        self.setupUi(self)
        self._setup_assignment_rows(station_ids, initial_records)


class DomainAndDischargeTableAssignmentDialog(
    QDialog,
    Ui_DomainAndDischargeTableAssignmentDialog,
    _OutletAssignmentRows,
):
    """Assign optional domains and discharge tables to supplied outlets."""

    include_domain = True

    def __init__(
        self,
        outlet_ids: list[str],
        parent=None,
        initial_records: Mapping[str, Mapping] | None = None,
    ):
        super().__init__(parent)
        self.setupUi(self)
        self._setup_assignment_rows(outlet_ids, initial_records)


__all__ = [
    "DischargeTableAssignmentDialog",
    "DomainAndDischargeTableAssignmentDialog",
    "OutletAssignment",
]

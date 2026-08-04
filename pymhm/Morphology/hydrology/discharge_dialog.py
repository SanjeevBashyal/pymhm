# -*- coding: utf-8 -*-
"""Discharge table assignment dialog for gauged outlets."""
from __future__ import annotations

try:
    from qgis.gui import QgsMapLayerComboBox
except Exception:  # pragma: no cover - fallback for generated UI imports
    from qgsmaplayercombobox import QgsMapLayerComboBox

from qgis.PyQt.QtWidgets import (
    QDialog,
    QLabel,
    QSizePolicy,
)

from ...qgis_compat import map_layer_filters
from ...pyui.ui_discharge_table_assignment_dialog import (
    Ui_DischargeTableAssignmentDialog,
)


class DischargeTableAssignmentDialog(QDialog, Ui_DischargeTableAssignmentDialog):
    """Dialog with one discharge table layer dropdown for each station ID."""

    def __init__(self, station_ids: list[str], parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self._combos: dict[str, QgsMapLayerComboBox] = {}

        for row, station_id in enumerate(station_ids):
            label = QLabel(f"STATION_ID: {station_id}", self.rows_widget)
            combo = QgsMapLayerComboBox(self.rows_widget)
            combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            combo.setFilters(map_layer_filters("VectorLayer"))
            if hasattr(combo, "setAllowEmptyLayer"):
                combo.setAllowEmptyLayer(False)
            self.rows_layout.addWidget(label, row, 0)
            self.rows_layout.addWidget(combo, row, 1)
            self._combos[station_id] = combo

    def selected_layers(self) -> dict[str, object]:
        """Return selected discharge table layers keyed by station ID."""
        return {
            station_id: combo.currentLayer()
            for station_id, combo in self._combos.items()
        }

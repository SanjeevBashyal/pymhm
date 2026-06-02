# -*- coding: utf-8 -*-
"""Discharge table assignment dialog for gauged outlets."""
from __future__ import annotations

try:
    from qgis.gui import QgsMapLayerComboBox
except Exception:  # pragma: no cover - fallback for generated UI imports
    from qgsmaplayercombobox import QgsMapLayerComboBox

from qgis.core import QgsMapLayerProxyModel
from qgis.PyQt.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QGridLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


class DischargeTableAssignmentDialog(QDialog):
    """Dialog with one discharge table layer dropdown for each station ID."""

    def __init__(self, station_ids: list[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Assign Discharge Tables")
        self.resize(560, 320)
        self._combos: dict[str, QgsMapLayerComboBox] = {}

        layout = QVBoxLayout(self)
        intro = QLabel("Select one discharge table layer for each STATION_ID.")
        intro.setWordWrap(True)
        layout.addWidget(intro)

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        rows_widget = QWidget(scroll_area)
        rows_layout = QGridLayout(rows_widget)
        rows_layout.setColumnStretch(1, 1)

        for row, station_id in enumerate(station_ids):
            label = QLabel(f"STATION_ID: {station_id}", rows_widget)
            combo = QgsMapLayerComboBox(rows_widget)
            combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            combo.setFilters(QgsMapLayerProxyModel.VectorLayer)
            if hasattr(combo, "setAllowEmptyLayer"):
                combo.setAllowEmptyLayer(False)
            rows_layout.addWidget(label, row, 0)
            rows_layout.addWidget(combo, row, 1)
            self._combos[station_id] = combo

        scroll_area.setWidget(rows_widget)
        layout.addWidget(scroll_area)

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected_layers(self) -> dict[str, object]:
        """Return selected discharge table layers keyed by station ID."""
        return {
            station_id: combo.currentLayer()
            for station_id, combo in self._combos.items()
        }

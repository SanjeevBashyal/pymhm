# -*- coding: utf-8 -*-
"""Elevation-band parameter dialog helpers."""
from __future__ import annotations

from ..common import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QVBoxLayout,
)
from ..core.naming import NamingAndRangeMixin


class ElevationBandDialogMixin(NamingAndRangeMixin):
    """Elevation-band parameter dialog helpers."""

    def _ask_elevation_band_width(
            self,
            min_elevation: float,
            max_elevation: float) -> float | None:
        """Show elevation range and ask for the elevation window width."""
        band_dialog = QDialog(self.dialog)
        band_dialog.setWindowTitle("Elevation Bands")

        layout = QVBoxLayout(band_dialog)
        form_layout = QFormLayout()

        form_layout.addRow(
            "Selected DEM minimum:",
            QLabel(f"{min_elevation:.2f}")
        )
        form_layout.addRow(
            "Selected DEM maximum:",
            QLabel(f"{max_elevation:.2f}")
        )

        elevation_range = max(max_elevation - min_elevation, 1.0)
        suggested_width = self._nice_step(elevation_range / 10.0)
        width_spin_box = QDoubleSpinBox()
        width_spin_box.setDecimals(3)
        width_spin_box.setMinimum(0.001)
        width_spin_box.setMaximum(max(elevation_range * 10.0, suggested_width))
        width_spin_box.setSingleStep(max(suggested_width / 10.0, 0.001))
        width_spin_box.setValue(suggested_width)
        form_layout.addRow("Elevation window width:", width_spin_box)

        layout.addLayout(form_layout)

        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(band_dialog.accept)
        button_box.rejected.connect(band_dialog.reject)
        layout.addWidget(button_box)

        if band_dialog.exec_() != QDialog.Accepted:
            return None

        return width_spin_box.value()

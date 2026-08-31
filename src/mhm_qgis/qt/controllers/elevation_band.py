# -*- coding: utf-8 -*-
"""Widget state for `elevation_band_dialog.ui`."""
from __future__ import annotations

try:
    from qgis.PyQt.QtWidgets import QDialog
except ImportError:
    from ...standalone import install

    install(force=True)
    from qgis.PyQt.QtWidgets import QDialog

from ...qt.ui.pyui.ui_elevation_band_dialog import Ui_ElevationBandDialog

def _ask_elevation_band_width(
        dialog,
        min_elevation: float,
        max_elevation: float) -> float | None:
    """Show elevation range and ask for the elevation window width."""
    band_dialog = QDialog(dialog.dialog)
    ui = Ui_ElevationBandDialog()
    ui.setupUi(band_dialog)
    ui.min_elevation_value.setText(f"{min_elevation:.2f}")
    ui.max_elevation_value.setText(f"{max_elevation:.2f}")

    elevation_range = max(max_elevation - min_elevation, 1.0)
    suggested_width = dialog._nice_step(elevation_range / 10.0)
    ui.width_spin_box.setMaximum(
        max(elevation_range * 10.0, suggested_width))
    ui.width_spin_box.setSingleStep(max(suggested_width / 10.0, 0.001))
    ui.width_spin_box.setValue(suggested_width)

    if band_dialog.exec_() != QDialog.Accepted:
        return None

    return ui.width_spin_box.value()

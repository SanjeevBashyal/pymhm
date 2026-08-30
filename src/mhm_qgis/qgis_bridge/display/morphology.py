# -*- coding: utf-8 -*-
"""Canvas half of `groupBox_displayMorphLayers`.

Which file to show is decided by the QGIS-free resolver in
:mod:`mhm_qgis.morphology_display`; this module only drives the canvas.
"""
from __future__ import annotations

from qgis.PyQt import QtCore
from qgis.PyQt.QtWidgets import QMessageBox

from ...morphology_display import land_cover_periods, resolve_display_output
from .canvas import show_raster


def update_date_control(dialog, index=None) -> None:
    """Enable the date selector only for historical land cover."""
    editor = dialog.dateTimeEdit_forHistoricalInputs
    combo = dialog.comboBox_morphVariableToDisplay
    if editor is None or combo is None:
        return
    key = combo.currentData()
    periods = land_cover_periods(dialog.project_folder) if dialog.project_folder else []
    enabled = key == "land_cover" and len(periods) > 1
    editor.setEnabled(enabled)
    if not enabled:
        return
    try:
        first = min(int(period["start_year"]) for period in periods)
        last = max(int(period["end_year"]) for period in periods)
    except (KeyError, TypeError, ValueError):
        editor.setEnabled(False)
        return
    editor.setMinimumDate(QtCore.QDate(first, 1, 1))
    editor.setMaximumDate(QtCore.QDate(last, 12, 31))
    editor.setDate(QtCore.QDate(first, 1, 1))


def display_selected_layer(dialog, checked=False) -> None:
    """Load the prepared output selected by the display combo box."""
    if not dialog.project_folder:
        QMessageBox.warning(dialog, "Display Layer", "Select a project folder first.")
        return
    key = dialog.comboBox_morphVariableToDisplay.currentData()
    if not key:
        QMessageBox.warning(dialog, "Display Layer", "Select a morphology layer.")
        return
    year = None
    if dialog.dateTimeEdit_forHistoricalInputs.isEnabled():
        year = dialog.dateTimeEdit_forHistoricalInputs.date().year()
    output = resolve_display_output(dialog.project_folder, key, year=year)
    if output is None:
        QMessageBox.warning(
            dialog,
            "Display Layer",
            "The selected morphology output has not been prepared yet.",
        )
        return
    show_raster(
        dialog,
        output.path,
        output.name,
        variable=output.variable,
        band=output.band if output.band and output.band > 1 else None,
        is_raster=output.is_raster,
    )


__all__ = ["display_selected_layer", "update_date_control"]

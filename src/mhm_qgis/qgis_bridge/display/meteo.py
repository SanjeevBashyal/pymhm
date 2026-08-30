# -*- coding: utf-8 -*-
"""Canvas half of `groupBox_meteoDisplay`.

Pick a forcing variable and a date, and show that timestep as a raster. The
slider scrubs the same period; both controls drive one shared step index.
"""
from __future__ import annotations

from qgis.PyQt import QtCore
from qgis.PyQt.QtWidgets import QMessageBox

from ...Meteorology.display import (
    VARIABLE_LABELS,
    available_meteo_variables,
    date_for_step,
    project_time_range,
    resolve_meteo_output,
    step_count,
    step_for_date,
)
from .canvas import show_raster

_DATE_FORMAT = "yyyy-MM-dd"


def _as_qdate(value) -> QtCore.QDate:
    return QtCore.QDate(value.year, value.month, value.day)


def _period(dialog):
    """Return the forcing period for the open project, if any."""
    if not dialog.project_folder:
        return None
    return project_time_range(dialog.project_folder)


def refresh(dialog) -> None:
    """Repopulate the variables and re-range the date and slider controls.

    Called after meteorology preparation and when a project is loaded, so the
    controls always describe what is actually on disk.
    """
    combo = dialog.comboBox_meteoVarDisplay
    slider = dialog.horizontalSlider_meteoTimeSelector
    editor = dialog.dateTimeEdit_meteoVarDisplay

    variables = (
        available_meteo_variables(dialog.project_folder)
        if dialog.project_folder
        else []
    )
    previous = combo.currentData()
    combo.blockSignals(True)
    combo.clear()
    for variable in variables:
        combo.addItem(VARIABLE_LABELS.get(variable, variable), variable)
    if previous in variables:
        combo.setCurrentIndex(variables.index(previous))
    combo.blockSignals(False)

    period = _period(dialog)
    enabled = bool(variables) and period is not None
    for widget in (combo, slider, editor, dialog.pushButton_meteoVarDisplay):
        widget.setEnabled(enabled)
    if not enabled:
        dialog.label_startTime.setText("")
        dialog.label_endTime.setText("")
        return

    first, last = period
    dialog.label_startTime.setText(first.isoformat())
    dialog.label_endTime.setText(last.isoformat())

    editor.blockSignals(True)
    editor.setDisplayFormat(_DATE_FORMAT)
    editor.setMinimumDate(_as_qdate(first))
    editor.setMaximumDate(_as_qdate(last))
    editor.setDate(_as_qdate(first))
    editor.blockSignals(False)

    slider.blockSignals(True)
    slider.setMinimum(0)
    slider.setMaximum(max(0, step_count(period) - 1))
    slider.setValue(0)
    slider.blockSignals(False)


def slider_moved(dialog, value: int) -> None:
    """Move the date editor to the step the slider selected."""
    period = _period(dialog)
    if period is None:
        return
    editor = dialog.dateTimeEdit_meteoVarDisplay
    editor.blockSignals(True)
    editor.setDate(_as_qdate(date_for_step(period, value)))
    editor.blockSignals(False)


def date_changed(dialog, value=None) -> None:
    """Move the slider to the step the date editor selected."""
    period = _period(dialog)
    if period is None:
        return
    chosen = dialog.dateTimeEdit_meteoVarDisplay.date().toPyDate()
    slider = dialog.horizontalSlider_meteoTimeSelector
    slider.blockSignals(True)
    slider.setValue(step_for_date(period, chosen))
    slider.blockSignals(False)


def display_selected_layer(dialog, checked=False) -> None:
    """Show the selected forcing variable at the selected date."""
    if not dialog.project_folder:
        QMessageBox.warning(dialog, "Display Meteorology", "Select a project folder first.")
        return
    variable = dialog.comboBox_meteoVarDisplay.currentData()
    if not variable:
        QMessageBox.warning(
            dialog,
            "Display Meteorology",
            "No prepared meteorology forcing was found for this project.",
        )
        return
    when = dialog.dateTimeEdit_meteoVarDisplay.date().toPyDate()
    output = resolve_meteo_output(dialog.project_folder, variable, when)
    if output is None:
        QMessageBox.warning(
            dialog,
            "Display Meteorology",
            f"{variable} has no data for {when.isoformat()}.",
        )
        return
    show_raster(
        dialog,
        output.path,
        output.name,
        variable=output.variable,
        band=output.band,
    )


__all__ = ["date_changed", "display_selected_layer", "refresh", "slider_moved"]

# -*- coding: utf-8 -*-
"""Widget wiring for `land_use_historical_input.ui` and `soil_multi_horizon_input.ui`.

Only the controls the form declares are wired here. The browse button of each
land-use or soil row is connected where that row is built, since it closes over
the row's own combo adapter.
"""
from __future__ import annotations


def bind_historical_land_use(dialog) -> None:
    """Wire the row-count control, the period table and the lookup selector."""
    dialog.pushButton_addLandUseInputWidgets.clicked.connect(
        lambda: dialog.set_layer_count(dialog.spinBox_nLandUseLayers.value())
    )
    dialog.tableWidget_landUseTimeInputs.cellChanged.connect(dialog._update_bounds)
    dialog.comboBox_lookupTableInput.currentIndexChanged.connect(dialog._lookup_changed)
    dialog.pushButton_browsecLookupTable.clicked.connect(dialog._browse_lookup)


def bind_multi_horizon_soil(dialog) -> None:
    """Wire the horizon-count control."""
    dialog.pushButton_addHorizonInputWidgets.clicked.connect(
        lambda: dialog.set_horizon_count(dialog.spinBox_nSoilHorizons.value())
    )


__all__ = ["bind_historical_land_use", "bind_multi_horizon_soil"]

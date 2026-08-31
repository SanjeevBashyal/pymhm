# -*- coding: utf-8 -*-
"""Widget wiring for `land_use_historical_input.ui`.

Each row's browse button is connected where the row is built, since it closes
over that row's combo adapter.
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


__all__ = ["bind_historical_land_use"]

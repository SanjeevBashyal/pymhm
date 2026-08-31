# -*- coding: utf-8 -*-
"""Widget wiring for `soil_multi_horizon_input.ui`."""
from __future__ import annotations


def bind_multi_horizon_soil(dialog) -> None:
    """Wire the horizon-count control."""
    dialog.pushButton_addHorizonInputWidgets.clicked.connect(
        lambda: dialog.set_horizon_count(dialog.spinBox_nSoilHorizons.value())
    )


__all__ = ["bind_multi_horizon_soil"]

# -*- coding: utf-8 -*-
"""Widget wiring for `mhm_ready_dialog.ui`."""
from __future__ import annotations

RASTER_FILTER = "Raster files (*.asc *.nc *.tif *.tiff);;All files (*)"
INPUT_FILTER = "Input layers (*.asc *.nc *.shp *.tif *.tiff);;All files (*)"
NETCDF_FILTER = "NetCDF (*.nc);;All files (*)"
CLASS_DEFINITION_FILTER = "Class definition (*.txt *.csv);;All files (*)"


def _browse():
    """Import lazily: `input_selection` imports this module at load time."""
    from ..dialogs.input_selection import _browse_into_combo

    return _browse_into_combo


def bind_mhm_ready(dialog, title: str) -> None:
    """Wire the layer and class-definition browse buttons."""
    browse = _browse()
    dialog.pushButton_browseInputLayer.clicked.connect(
        lambda: browse(
            dialog,
            dialog.comboBox_inputLayer,
            f"Select mHM-ready {title}",
            RASTER_FILTER,
        )
    )
    dialog.pushButton_browseClassDefinition.clicked.connect(
        lambda: browse(
            dialog,
            dialog.comboBox_classDefinitionInput,
            "Select class definition",
            CLASS_DEFINITION_FILTER,
        )
    )


__all__ = ["bind_mhm_ready"]

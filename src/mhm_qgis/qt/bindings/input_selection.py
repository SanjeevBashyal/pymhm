# -*- coding: utf-8 -*-
"""Widget wiring for the input-selection forms.

Covers `single_layer_input_with_lookup_config_dialog.ui`, `mhm_ready_dialog.ui`
and `lai_netcdf_input_dialog.ui`. Several browse buttons need the dialog title
chosen at construction, so it is passed in rather than rediscovered here.
"""
from __future__ import annotations

RASTER_FILTER = "Raster files (*.asc *.nc *.tif *.tiff);;All files (*)"
INPUT_FILTER = "Input layers (*.asc *.nc *.shp *.tif *.tiff);;All files (*)"
NETCDF_FILTER = "NetCDF (*.nc);;All files (*)"
CLASS_DEFINITION_FILTER = "Class definition (*.txt *.csv);;All files (*)"


def _browse():
    """Import lazily: `input_selection` imports this module at load time."""
    from ...input_selection import _browse_into_combo

    return _browse_into_combo


def bind_lookup_config(dialog) -> None:
    """Wire the lookup table and its two field selectors."""
    dialog.lookup_table_combo.currentIndexChanged.connect(dialog._refresh_fields)
    dialog.mapping_field_combo.currentIndexChanged.connect(dialog._update_ok)
    dialog.class_field_combo.currentIndexChanged.connect(dialog._update_ok)


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


def bind_single_layer_lookup(dialog, title: str) -> None:
    """Wire the input layer, lookup table and their refresh."""
    browse = _browse()
    dialog.lookup_table_combo.currentIndexChanged.connect(dialog._refresh_fields)
    dialog.pushButton_browseInputLayer.clicked.connect(
        lambda: browse(
            dialog,
            dialog.comboBox_inputLayer,
            f"Select {title} input",
            INPUT_FILTER,
        )
    )
    dialog.pushButton_browseLookupTable.clicked.connect(dialog._browse_lookup)


def bind_lai_netcdf(dialog) -> None:
    """Wire the LAI NetCDF browse button."""
    browse = _browse()
    dialog.pushButton_browseLAINetCDFInput.clicked.connect(
        lambda: browse(
            dialog,
            dialog.comboBox_laiNetCDFInput,
            "Select LAI NetCDF",
            NETCDF_FILTER,
        )
    )


__all__ = [
    "bind_lai_netcdf",
    "bind_lookup_config",
    "bind_mhm_ready",
    "bind_single_layer_lookup",
]

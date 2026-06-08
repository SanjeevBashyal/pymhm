# -*- coding: utf-8 -*-
"""Compatibility helpers for QGIS API changes used by the plugin UI."""

from qgis.core import Qgis, QgsMapLayerProxyModel


def map_layer_filters(*filter_names):
    """Return map layer filter flags using the modern Qgis.LayerFilter enum."""
    layer_filter_enum = getattr(Qgis, "LayerFilter", None)
    if layer_filter_enum is not None:
        filters = None
        for filter_name in filter_names:
            value = getattr(layer_filter_enum, filter_name)
            filters = value if filters is None else filters | value
        if filters is not None:
            return filters

    legacy_filters = None
    for filter_name in filter_names:
        value = getattr(QgsMapLayerProxyModel, filter_name)
        legacy_filters = value if legacy_filters is None else legacy_filters | value
    return legacy_filters

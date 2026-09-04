# -*- coding: utf-8 -*-
"""QGIS layer handling: build a layer, register it, find it, ask its CRS.

The one place in the plugin that knows the `QgsRasterLayer` / `QgsVectorLayer` /
`QgsProject` vocabulary. Everything above it -- the Qt dialogs, the display
helpers, the Morphology processors -- goes through this API instead of touching
the QGIS layer classes directly, so the version shims and the validity checks
live once.

QGIS is imported inside the functions, never at module scope, so `qgis_bridge`
stays importable in tests and under the standalone shim.
"""
from .compat import create_vector_file_writer, map_layer_filters, qgs_field
from .crs import crs_of, transform_between, transform_to_raster
from .gauges import gauge_coordinates
from .loader import load, open_layer, raster_source, source_uri
from .materialize import (
    local_layer_source,
    materialize_raster_layer,
    materialize_vector_layer,
    remove_vector_dataset,
)
from .project import (
    add,
    find_by_name,
    find_by_source,
    normalized_source,
    remove_under,
)

__all__ = [
    "add",
    "create_vector_file_writer",
    "crs_of",
    "transform_between",
    "find_by_name",
    "find_by_source",
    "gauge_coordinates",
    "load",
    "local_layer_source",
    "map_layer_filters",
    "normalized_source",
    "materialize_vector_layer",
    "materialize_raster_layer",
    "open_layer",
    "qgs_field",
    "raster_source",
    "remove_under",
    "remove_vector_dataset",
    "source_uri",
    "transform_to_raster",
]

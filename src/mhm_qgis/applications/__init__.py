# -*- coding: utf-8 -*-
"""Handlers for the applications the plugin depends on but does not own.

One module per external package. Everything the plugin asks of mhm-tools or
nml-tools goes through these, so there is a single place that names them and a
single place to look when one of them changes.

Their imports stay inside the functions: importing ``mhm_tools`` costs roughly
850 ms and pulls in geopandas, pandas, rasterio, scipy and xarray, which must
not be paid when QGIS loads the plugin.
"""

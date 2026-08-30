# -*- coding: utf-8 -*-
"""Running the plugin outside QGIS: the Qt-backed shim and its entry point."""
from .qgis_shim import install, is_active

__all__ = ["install", "is_active"]

# -*- coding: utf-8 -*-
"""Meteorology processing package for pymhm."""

__all__ = ["MeteorologyProcessor"]


def __getattr__(name):
    """Keep pure meteorology helpers importable without a QGIS runtime."""
    if name == "MeteorologyProcessor":
        from .processor import MeteorologyProcessor

        return MeteorologyProcessor
    raise AttributeError(name)

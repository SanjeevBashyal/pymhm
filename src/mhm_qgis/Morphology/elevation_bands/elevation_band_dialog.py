# -*- coding: utf-8 -*-
"""Elevation-band parameter dialog helpers."""
from __future__ import annotations

from ..common import QDialog
from ..core.naming import NamingAndRangeMixin
from ...qt.controllers import elevation_band as elevation_band_controller
from ...qt.ui.pyui.ui_elevation_band_dialog import Ui_ElevationBandDialog


class ElevationBandDialogMixin(NamingAndRangeMixin):
    """Elevation-band parameter dialog helpers."""

    def _ask_elevation_band_width(self, *args, **kwargs):
        return elevation_band_controller._ask_elevation_band_width(self, *args, **kwargs)

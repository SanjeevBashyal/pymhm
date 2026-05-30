# -*- coding: utf-8 -*-
"""Compatibility aggregate for hydrological morphology mixins."""
from .terrain import TerrainAnalysisMixin
from .flow import FlowAnalysisMixin
from .gauge import GaugePositionMixin


class HydrologyMixin(
    TerrainAnalysisMixin,
    FlowAnalysisMixin,
    GaugePositionMixin,
):
    """Compatibility aggregate for hydrological morphology mixins."""

    pass

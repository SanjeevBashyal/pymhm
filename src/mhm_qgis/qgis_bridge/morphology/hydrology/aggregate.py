# -*- coding: utf-8 -*-
"""Compatibility aggregate for hydrological morphology mixins."""

from .flow import FlowAnalysisMixin
from .gauge import GaugePositionMixin
from .terrain import TerrainAnalysisMixin


class HydrologyMixin(
    GaugePositionMixin,
    TerrainAnalysisMixin,
    FlowAnalysisMixin,
):
    """Compatibility aggregate for hydrological morphology mixins."""

    pass

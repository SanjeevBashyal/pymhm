# -*- coding: utf-8 -*-
"""Compatibility aggregate for hydrological morphology mixins."""
from .terrain import TerrainAnalysisMixin
from .flow import FlowAnalysisMixin
from .gauge import GaugePositionMixin
from .outlets import OutletCountMixin
from .discharge import DischargeAssignmentMixin


class HydrologyMixin(
    TerrainAnalysisMixin,
    FlowAnalysisMixin,
    GaugePositionMixin,
    OutletCountMixin,
    DischargeAssignmentMixin,
):
    """Compatibility aggregate for hydrological morphology mixins."""

    pass

# -*- coding: utf-8 -*-
"""Compatibility aggregate for hydrological morphology mixins."""

from .discharge import DischargeAssignmentMixin
from .flow import FlowAnalysisMixin
from .gauge import GaugePositionMixin
from .outlets import OutletCountMixin
from .terrain import TerrainAnalysisMixin


class HydrologyMixin(
    GaugePositionMixin,
    TerrainAnalysisMixin,
    FlowAnalysisMixin,
    OutletCountMixin,
    DischargeAssignmentMixin,
):
    """Compatibility aggregate for hydrological morphology mixins."""

    pass

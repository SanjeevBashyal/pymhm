# -*- coding: utf-8 -*-
"""Compatibility aggregate for watershed morphology mixins."""
from .dem_fill import DemFillMixin
from .channel_network import ChannelNetworkMixin
from .pour_points import PourPointMixin
from .watershed_delineation import WatershedDelineationMixin


class WatershedMixin(
    WatershedDelineationMixin,
    PourPointMixin,
    ChannelNetworkMixin,
    DemFillMixin,
):
    """Compatibility aggregate for watershed morphology mixins."""

    pass

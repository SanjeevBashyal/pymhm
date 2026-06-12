# -*- coding: utf-8 -*-
"""Compatibility aggregate for layer-processing morphology mixins."""
from .land_cover import LandCoverProcessingMixin
from .soil import SoilProcessingMixin
from .geology import GeologyProcessingMixin
from .masking import MaskingMixin


class LayerProcessingMixin(
    SoilProcessingMixin,
    GeologyProcessingMixin,
    LandCoverProcessingMixin,
    MaskingMixin,
):
    """Compatibility aggregate for layer-processing morphology mixins."""

    pass

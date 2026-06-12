# -*- coding: utf-8 -*-
"""Compatibility aggregate for layer-processing morphology mixins."""
from .land_cover import LandCoverProcessingMixin
from .soil import SoilProcessingMixin
from .geology import GeologyProcessingMixin
from .lai import LaiProcessingMixin
from .masking import MaskingMixin


class LayerProcessingMixin(
    SoilProcessingMixin,
    GeologyProcessingMixin,
    LandCoverProcessingMixin,
    LaiProcessingMixin,
    MaskingMixin,
):
    """Compatibility aggregate for layer-processing morphology mixins."""

    pass

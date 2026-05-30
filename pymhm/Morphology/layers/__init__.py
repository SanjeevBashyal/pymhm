# -*- coding: utf-8 -*-
"""Compatibility aggregate for layer-processing morphology mixins."""
from .layer_preparation import LayerPreparationMixin
from .land_cover_lookup import LandCoverLookupMixin
from .land_cover import LandCoverProcessingMixin
from .soil import SoilProcessingMixin
from .geology import GeologyProcessingMixin
from .masking import MaskingMixin


class LayerProcessingMixin(
    LayerPreparationMixin,
    LandCoverLookupMixin,
    LandCoverProcessingMixin,
    SoilProcessingMixin,
    GeologyProcessingMixin,
    MaskingMixin,
):
    """Compatibility aggregate for layer-processing morphology mixins."""

    pass

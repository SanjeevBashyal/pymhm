# -*- coding: utf-8 -*-
"""Compatibility aggregate for layer-processing morphology mixins."""
from .categorical import CategoricalProcessingMixin
from .lai import LaiProcessingMixin
from .masking import MaskingMixin


class LayerProcessingMixin(
    CategoricalProcessingMixin,
    LaiProcessingMixin,
    MaskingMixin,
):
    """Compatibility aggregate for layer-processing morphology mixins."""

    pass

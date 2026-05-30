# -*- coding: utf-8 -*-
"""Composed morphology/geometry processor for pymhm."""
from .common import DialogUtils
from .core import CoreMixin
from .watershed import WatershedMixin
from .elevation_bands import ElevationBandsMixin
from .hydrology import HydrologyMixin
from .layers import LayerProcessingMixin
from .latlon import LatLonMixin
from .classification_writers import ClassificationWritersMixin
from .orchestration import OrchestrationMixin


class MorphologyProcessor(
    CoreMixin,
    WatershedMixin,
    ElevationBandsMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
    ClassificationWritersMixin,
    OrchestrationMixin,
    DialogUtils,
):
    """Handles all morphology/geometry processing functionality."""

    pass

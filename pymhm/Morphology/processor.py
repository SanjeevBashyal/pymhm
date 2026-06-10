# -*- coding: utf-8 -*-
"""Composed morphology/geometry processor for pymhm."""
from .common import DialogUtils
from .core.project_state import ProjectStateMixin
from .watershed import WatershedMixin
from .elevation_bands import ElevationBandsMixin
from .hydrology.aggregate import HydrologyMixin
from .layers import LayerProcessingMixin
from .latlon import LatLonMixin
from .classification_writers import ClassificationWritersMixin
from .orchestration.execute_all import ExecuteAllMixin
from .orchestration.reset_geometry import ResetGeometryMixin


class MorphologyProcessor(
    ElevationBandsMixin,
    ProjectStateMixin,
    ExecuteAllMixin,
    WatershedMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
    ClassificationWritersMixin,
    ResetGeometryMixin,
    DialogUtils,
):
    """Handles all morphology/geometry processing functionality."""

    pass

# -*- coding: utf-8 -*-
"""Composed morphology/geometry processor for pymhm."""
from .common import DialogUtils
from .core.project_state import ProjectStateMixin
from .elevation_bands import ElevationBandsMixin
from .hydrology.aggregate import HydrologyMixin
from .latlon import LatLonMixin
from .layers import LayerProcessingMixin
from .orchestration.execute_all import ExecuteAllMixin
from .orchestration.reset_geometry import ResetGeometryMixin
from .watershed import WatershedMixin


class MorphologyProcessor(
    ElevationBandsMixin,
    ProjectStateMixin,
    ExecuteAllMixin,
    WatershedMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
    ResetGeometryMixin,
    DialogUtils,
):
    """Handles all morphology/geometry processing functionality."""

    pass

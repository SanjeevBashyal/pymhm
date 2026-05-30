# -*- coding: utf-8 -*-
"""Compatibility aggregate for core morphology infrastructure mixins."""
from .base import BaseProcessingMixin
from .processing_state import ProcessingStateMixin
from .project_state import ProjectStateMixin
from .dem_inputs import DemInputMixin
from .dependencies import PythonDependencyMixin
from .raster_io import RasterIOMixin
from .vector_io import VectorIOMixin
from .naming import NamingAndRangeMixin
from .predecessors import PredecessorMixin


class CoreMixin(
    BaseProcessingMixin,
    ProcessingStateMixin,
    ProjectStateMixin,
    DemInputMixin,
    PythonDependencyMixin,
    RasterIOMixin,
    VectorIOMixin,
    NamingAndRangeMixin,
    PredecessorMixin,
):
    """Compatibility aggregate for core morphology infrastructure mixins."""

    pass

# -*- coding: utf-8 -*-
"""Compatibility aggregate for core morphology infrastructure mixins."""
from .project_state import ProjectStateMixin
from .dem_inputs import DemInputMixin
from .raster_io import RasterIOMixin
from .vector_io import VectorIOMixin
from .predecessors import PredecessorMixin


class CoreMixin(
    ProjectStateMixin,
    DemInputMixin,
    RasterIOMixin,
    VectorIOMixin,
    PredecessorMixin,
):
    """Compatibility aggregate for core morphology infrastructure mixins."""

    pass

# -*- coding: utf-8 -*-
"""Compatibility aggregate for pour-point snapping workflows."""
from .pour_point_workflow import PourPointWorkflowMixin
from .network_snapper import NetworkSnapperMixin


class PourPointMixin(
    PourPointWorkflowMixin,
):
    """Compatibility aggregate for pour-point snapping workflows."""

    pass

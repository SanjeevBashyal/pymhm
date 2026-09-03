# -*- coding: utf-8 -*-
"""Compatibility aggregate for reset and execute-all morphology workflows."""
from .reset_geometry import ResetGeometryMixin
from .execute_all import ExecuteAllMixin


class OrchestrationMixin(
    ExecuteAllMixin,
    ResetGeometryMixin,
):
    """Compatibility aggregate for reset and execute-all morphology workflows."""

    pass

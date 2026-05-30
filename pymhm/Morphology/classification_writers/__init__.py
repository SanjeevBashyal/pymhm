# -*- coding: utf-8 -*-
"""Compatibility aggregate for morphology classdefinition writers."""
from .geology_classification import GeologyClassificationWriterMixin
from .soil_classification import SoilClassificationWriterMixin


class ClassificationWritersMixin(
    GeologyClassificationWriterMixin,
    SoilClassificationWriterMixin,
):
    """Compatibility aggregate for morphology classdefinition writers."""

    pass

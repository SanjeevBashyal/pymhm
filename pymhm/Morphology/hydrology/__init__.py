# -*- coding: utf-8 -*-
"""Hydrological morphology mixins."""

__all__ = ["HydrologyMixin"]


def __getattr__(name):
    """Load the aggregate mixin lazily to avoid submodule import cycles."""
    if name == "HydrologyMixin":
        from .aggregate import HydrologyMixin
        return HydrologyMixin
    raise AttributeError(name)

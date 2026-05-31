# -*- coding: utf-8 -*-
"""Meteorology workflow exceptions."""
from __future__ import annotations


class MeteorologyInputError(ValueError):
    """Raised when a meteorology run cannot start because inputs are missing."""

    def __init__(self, title: str, message: str, severity: str = "warning"):
        super().__init__(message)
        self.title = title
        self.message = message
        self.severity = severity


class MeteorologyToolImportError(RuntimeError):
    """Raised when the ERA5-Land processing backend cannot be imported."""

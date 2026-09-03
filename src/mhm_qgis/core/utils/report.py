# -*- coding: utf-8 -*-
"""Severity levels for the `report=` callback contract.

Defined here, in a leaf with no dependencies, so `core` can name a severity
without importing the Qt layer that renders it. `qt/reporting` imports these;
nothing goes the other way.
"""
from __future__ import annotations

WARNING = "warning"
CRITICAL = "critical"
INFORMATION = "information"

__all__ = ["CRITICAL", "INFORMATION", "WARNING"]

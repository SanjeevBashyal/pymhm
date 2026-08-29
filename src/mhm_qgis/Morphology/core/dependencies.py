# -*- coding: utf-8 -*-
"""Lazy loading and validation of optional Python morphology dependencies."""
from ..common import (
    QMessageBox,
    processing,
)


class PythonDependencyMixin:
    """Lazy loading and validation of optional Python morphology dependencies."""

    def _get_python_morphology_deps(self):
        """
        Load the Python morphology stack used in the workshop notebook.

        Imports are lazy so the plugin can still open even if a user's QGIS
        Python environment has not installed the optional morphology packages.
        """
        if hasattr(self, "_python_morphology_deps"):
            return self._python_morphology_deps

        missing = []
        deps = {}

        try:
            import numpy as np
            deps["np"] = np
        except ImportError:
            missing.append("numpy")

        try:
            import pyflwdir as pfd
            deps["pfd"] = pfd
        except ImportError:
            missing.append("pyflwdir")

        try:
            from affine import Affine
            deps["Affine"] = Affine
        except ImportError:
            missing.append("affine")

        try:
            from osgeo import gdal, osr
            deps["gdal"] = gdal
            deps["osr"] = osr
        except ImportError:
            missing.append("GDAL Python bindings (osgeo)")

        if missing:
            message = (
                "Python morphology processing requires these packages in the "
                f"QGIS Python environment: {', '.join(missing)}."
            )
            self.log_message(f"ERROR: {message}")
            QMessageBox.critical(
                self.dialog, "Missing Python Dependencies", message)
            return None

        self._python_morphology_deps = deps
        return deps

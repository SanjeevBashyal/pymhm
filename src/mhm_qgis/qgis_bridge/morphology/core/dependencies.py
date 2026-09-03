# -*- coding: utf-8 -*-
"""Availability of the optional Python morphology stack, reported to the user."""
from ....applications import pyflwdir_handler
from ....qt import reporting


class PythonDependencyMixin:
    """Availability of the optional Python morphology stack."""

    def _get_python_morphology_deps(self):
        """Return the morphology stack, or None after telling the user what is missing.

        The import list lives in `applications/pyflwdir_handler`; this only adds
        the caching and the user-facing message.
        """
        if hasattr(self, "_python_morphology_deps"):
            return self._python_morphology_deps

        try:
            deps = pyflwdir_handler.dependencies()
        except pyflwdir_handler.MissingDependencies as error:
            reporting.dialog_reporter(
                self.dialog, log=self.log_message
            )(reporting.CRITICAL, "Missing Python Dependencies", str(error))
            return None

        self._python_morphology_deps = deps
        return deps

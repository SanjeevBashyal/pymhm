# -*- coding: utf-8 -*-
"""The worker that runs the QGIS-free meteorology phase off the GUI thread."""
from __future__ import annotations

import traceback

try:
    from qgis.PyQt import QtCore
except ImportError:
    from ...standalone import install

    install(force=True)
    from qgis.PyQt import QtCore


class MeteorologyWorkflowWorker(QtCore.QObject):
    """Run only the QGIS-free meteorology phase in a worker thread."""

    log_message = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(str, bool, str)

    def __init__(
            self,
            workflow_key,
            workflow_label,
            meteorology_processor,
            meteorology_run):
        super().__init__()
        self.workflow_key = workflow_key
        self.workflow_label = workflow_label
        self.meteorology_processor = meteorology_processor
        self.meteorology_run = meteorology_run
        self._original_meteo_log_message = None

    @QtCore.pyqtSlot()
    def run(self):
        """Execute meteorology and hand morphology back to the GUI thread."""
        self._original_meteo_log_message = self.meteorology_processor.log_message
        self.meteorology_processor.log_message = self.log_message.emit
        try:
            ok = self.meteorology_processor.process_meteo_forcing(
                self.meteorology_run,
                show_dialog=False,
            )
            self.finished.emit(self.workflow_key, ok, "")
        except Exception as exc:
            self.log_message.emit(
                f"\nERROR: {self.workflow_label} worker failed: {exc}"
            )
            self.log_message.emit(f"Traceback: {traceback.format_exc()}")
            self.finished.emit(self.workflow_key, False, str(exc))
        finally:
            self.meteorology_processor.log_message = (
                self._original_meteo_log_message
            )

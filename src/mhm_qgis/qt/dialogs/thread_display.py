"""Task monitor dialog for the plugin worker queue."""
from __future__ import annotations

try:
    from qgis.PyQt import QtWidgets
except ImportError:
    from ...standalone import install

    install(force=True)
    from qgis.PyQt import QtWidgets

from ...qt.controllers import thread_display as thread_display_controller
from ...qt.bindings.thread_display import bind as bind_thread_display
from ...qt.ui.pyui.ui_thread_display_dialog import Ui_ThreadDisplayDialog


class ThreadDisplayDialog(QtWidgets.QDialog, Ui_ThreadDisplayDialog):
    def __init__(self, coordinator, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.coordinator = coordinator
        self.spinBox_threadCount.setValue(coordinator.max_threads)
        bind_thread_display(self, coordinator)
        coordinator.changed.connect(self.refresh)
        self.refresh()

    def refresh(self, *args, **kwargs):
        return thread_display_controller.refresh(self, *args, **kwargs)

    def _show_log(self, *args, **kwargs):
        return thread_display_controller._show_log(self, *args, **kwargs)


__all__ = ["ThreadDisplayDialog"]

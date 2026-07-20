"""Regression tests for UI forwarding of captured tool messages."""

import io
import logging
import sys

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from pymhm.mhm_tools_to_integrate.logging import capture_messages  # noqa: E402
from pymhm.utils import DialogUtils  # noqa: E402


class _LogText:
    def __init__(self):
        self.messages = []

    def append(self, message):
        self.messages.append(message)


class _Dialog(DialogUtils):
    def __init__(self):
        self.LogText = _LogText()


def test_captured_info_is_forwarded_once(monkeypatch):
    """A UI console echo must not feed back into redirected stdout."""
    console = io.StringIO()
    monkeypatch.setattr(sys, "__stdout__", console)
    dialog = _Dialog()

    with capture_messages(dialog.log_message):
        logging.getLogger("mhm_tools.pre.test").info("definition written")

    assert dialog.LogText.messages == ["INFO: definition written"]
    assert console.getvalue() == "INFO: definition written\n"

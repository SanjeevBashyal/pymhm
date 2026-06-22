# -*- coding: utf-8 -*-
"""Logging and stdout/stderr capture helpers for UI and CLI callers."""
from __future__ import annotations

import contextlib
import io
import logging
from collections.abc import Callable, Generator


LogCallback = Callable[[str], None]


class _StreamToCallback(io.TextIOBase):
    """File-like stream that forwards completed lines to a callback."""

    def __init__(self, callback: LogCallback, prefix: str = "") -> None:
        self._callback = callback
        self._prefix = prefix
        self._buffer = ""

    def writable(self) -> bool:
        return True

    def write(self, text: str) -> int:
        if not text:
            return 0
        self._buffer += str(text)
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            self._emit(line)
        return len(text)

    def flush(self) -> None:
        if self._buffer:
            self._emit(self._buffer)
            self._buffer = ""

    def _emit(self, line: str) -> None:
        line = line.rstrip("\r")
        if line:
            self._callback(f"{self._prefix}{line}")


class _CallbackLoggingHandler(logging.Handler):
    """Logging handler that forwards records to a callback."""

    def __init__(self, callback: LogCallback) -> None:
        super().__init__()
        self._callback = callback
        self.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._callback(self.format(record))
        except Exception:
            self.handleError(record)


@contextlib.contextmanager
def capture_messages(
    callback: LogCallback | None,
    logger_names: tuple[str, ...] = ("mhm_tools_to_integrate", "mhm_tools"),
) -> Generator[None, None, None]:
    """Capture stdout, stderr, and selected Python loggers into callback."""
    if callback is None:
        yield
        return

    stdout = _StreamToCallback(callback)
    stderr = _StreamToCallback(callback, prefix="ERROR: ")
    handlers: list[tuple[logging.Logger, logging.Handler, int, bool]] = []
    try:
        for logger_name in logger_names:
            logger = logging.getLogger(logger_name)
            handler = _CallbackLoggingHandler(callback)
            handlers.append((logger, handler, logger.level, logger.propagate))
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
            logger.propagate = False

        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            yield
    finally:
        stdout.flush()
        stderr.flush()
        for logger, handler, old_level, old_propagate in handlers:
            logger.removeHandler(handler)
            logger.setLevel(old_level)
            logger.propagate = old_propagate

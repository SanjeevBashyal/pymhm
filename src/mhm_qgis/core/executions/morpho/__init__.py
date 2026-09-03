# -*- coding: utf-8 -*-
"""What `pushButton_executeAllMorphology` runs, and in what order.

The sequence is data, and `run_stages` walks it. Neither knows about QGIS,
`QgsTask` or the dialog: the caller supplies a `start` function that knows how
to launch one stage and call back when it finishes. That keeps the answer to
"which stages, in what order, with what options" in one readable place instead
of buried in a list of closures inside the task bridge.

There used to be a second, synchronous implementation of this sequence that no
button reached. Keeping two in step was a standing trap -- a stage added to only
one silently never ran, which is how the LAI stage was once missed -- so the
unreachable twin is gone and this is the only definition.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Stage:
    """One step of Execute All: which starter to call, and with what."""

    starter: str
    args: tuple = ()
    options: dict = field(default_factory=dict)
    note: str = ""


#: Execute All, in order. `facc` and `fdir` come out of the derivatives stage,
#: so the hydrology pass only adds the upstream-area raster and the channel
#: network -- hence `accumulation=False, direction=False`.
EXECUTE_ALL_STAGES = (
    Stage("dem_derivatives", options={"key": "execute-all-dem"}),
    Stage("categorical", args=("lc",),
          options={"key": "execute-all-lc", "reuse_existing": True}),
    Stage("categorical", args=("soil",),
          options={"key": "execute-all-soil", "reuse_existing": True}),
    Stage("categorical", args=("geology",),
          options={"key": "execute-all-geology"}),
    Stage("lai"),
    Stage("hydrology",
          options={"key": "execute-all-hydrology", "load": "",
                   "accumulation": False, "direction": False},
          note="facc and fdir already came from the derivatives stage"),
)


def run_stages(stages, start, *, on_done, on_fail):
    """Walk `stages` one at a time, each continuing the next when it completes.

    `start(stage, advance)` launches one stage and returns falsy if it could not
    be started. Stages run asynchronously, so this returns as soon as the first
    one is under way; `on_done` fires after the last.
    """
    stages = tuple(stages)

    def advance(index=0):
        if index >= len(stages):
            on_done()
            return
        if not start(stages[index], lambda: advance(index + 1)):
            on_fail("Could not start the next preprocessing task.")

    advance()
    return True


__all__ = ["EXECUTE_ALL_STAGES", "Stage", "run_stages"]

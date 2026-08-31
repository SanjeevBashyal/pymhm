# -*- coding: utf-8 -*-
"""Widget state for each form: read, populate, enable, validate.

One module per form in `qt/ui/forms`, named after it, mirroring `qt/bindings`.
Each exposes functions taking the dialog as their first argument, so the
behaviour is testable without constructing the dialog's whole collaborator
graph.

The split against its neighbours:

- `qt/bindings` connects a signal to a slot and nothing else.
- `qt/controllers` (here) owns what a widget shows and whether it is enabled.
- `qgis_bridge` drives the QGIS application itself -- canvas, layers, project.
- `core` computes without knowing Qt or QGIS exist.

A dialog keeps a one-line method for any name another package calls, because
the dialog is the object processors and bridges hold a reference to.
"""

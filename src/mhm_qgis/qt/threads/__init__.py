# -*- coding: utf-8 -*-
"""Work that runs off the GUI thread.

`QThread` plus a `QObject` worker, for the long meteorology phase the dialog
must not block on. Distinct from `TaskCoordinator`/`QgsTask`, which QGIS
schedules for path-only raster work; this is the plain-Qt route used where the
plugin owns the thread itself.

Widgets, `QgsProject` and `QgsMapLayer` stay on the main thread. A worker gets
copied primitives and paths, and reports back through signals.
"""

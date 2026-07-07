# -*- coding: utf-8 -*-
"""QGIS plugin entry point for pymhm."""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load pymhm class from file pymhm.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """

    parent = iface.mainWindow() if iface is not None else None
    try:
        from .dependency_bootstrap import ensure_qgis_runtime_dependencies

        ensure_qgis_runtime_dependencies(parent=parent, prompt=True)
    except Exception as exc:
        try:
            from qgis.PyQt.QtWidgets import QMessageBox

            QMessageBox.warning(
                parent,
                "PymHM Dependency Check",
                f"PymHM dependency check could not run.\n\n{exc}",
            )
        except Exception:
            print(f"PymHM dependency check could not run: {exc}")

    from .pymhm import pymhm

    return pymhm(iface)

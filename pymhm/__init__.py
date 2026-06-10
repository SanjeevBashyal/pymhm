# -*- coding: utf-8 -*-
"""QGIS plugin entry point for pymhm."""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load pymhm class from file pymhm.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """

    from .pymhm import pymhm

    return pymhm(iface)

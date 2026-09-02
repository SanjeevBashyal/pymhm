# -*- coding: utf-8 -*-
"""Shared imports for mhm_qgis morphology processing modules."""
import os
import math
import csv
import json

from qgis.PyQt.QtWidgets import (
    QMessageBox,
    QFileDialog,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QVBoxLayout,
)
from qgis.core import (
    QgsVectorLayer,
    QgsRasterLayer,
    QgsApplication,
    QgsFeature,
    QgsGeometry,
    QgsSpatialIndex,
    QgsFields,
    QgsWkbTypes,
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsRectangle,
    QgsPointXY,
)
from qgis.PyQt.QtCore import NULL
import processing

from ..utils import DialogUtils
from ..core.handlers.store.paths import geometry_folder, geometry_folder as project_geometry_folder, morph_staging_folder
from ..core.handlers.store.layout import morph_folder

# -*- coding: utf-8 -*-
"""Shared imports for pymhm morphology processing modules."""
import os
import math
import csv
import json
from datetime import datetime

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
    QgsVectorFileWriter,
    QgsField,
    QgsFields,
    QgsWkbTypes,
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsRectangle,
    QgsPointXY,
)
from qgis.PyQt.QtCore import QVariant, NULL
import processing

from ..utils import DialogUtils
from ..project_layout import (
    geometry_folder,
    geometry_folder as project_geometry_folder,
    morph_folder,
)

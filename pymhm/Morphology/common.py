# -*- coding: utf-8 -*-
"""Shared imports for pymhm morphology processing modules."""
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
from qgis.PyQt.QtCore import NULL, QMetaType
try:
    from qgis.PyQt.QtCore import QVariant
except ImportError:
    QVariant = None
import processing

from ..utils import DialogUtils
from ..project_layout import (
    geometry_folder,
    geometry_folder as project_geometry_folder,
    morph_folder,
)


def _qmeta_type(type_name):
    """Return a QMetaType enum value across Qt/QGIS bindings."""
    enum_class = getattr(QMetaType, "Type", QMetaType)
    return getattr(enum_class, type_name)


def qgs_field(name, field_type):
    """Create a QgsField across old QVariant and newer QMetaType bindings."""
    qmeta_name_by_qvariant_name = {
        "String": "QString",
        "Int": "Int",
        "Double": "Double",
    }
    qvariant_type_by_name = {}
    if QVariant is not None:
        qvariant_type_by_name = {
            "String": QVariant.String,
            "Int": QVariant.Int,
            "Double": QVariant.Double,
        }
    qmeta_name = qmeta_name_by_qvariant_name.get(field_type, field_type)
    candidates = [_qmeta_type(qmeta_name)]
    try:
        candidates.append(int(candidates[0]))
    except (TypeError, ValueError):
        pass
    if field_type in qvariant_type_by_name:
        candidates.append(qvariant_type_by_name[field_type])

    last_error = None
    for candidate in candidates:
        try:
            return QgsField(name, candidate)
        except TypeError as exc:
            last_error = exc

    if last_error is not None:
        raise last_error
    return QgsField(name)


def create_vector_file_writer(path, fields, geometry_type, crs,
                              driver_name="ESRI Shapefile"):
    """Create a QgsVectorFileWriter using the non-deprecated factory when present."""
    if hasattr(QgsVectorFileWriter, "create"):
        options = QgsVectorFileWriter.SaveVectorOptions()
        options.driverName = driver_name
        options.fileEncoding = "UTF-8"
        try:
            writer = QgsVectorFileWriter.create(
                path,
                fields,
                geometry_type,
                crs,
                QgsProject.instance().transformContext(),
                options,
            )
            if isinstance(writer, tuple):
                return writer[0]
            return writer
        except TypeError:
            pass

    return QgsVectorFileWriter(
        path,
        "UTF-8",
        fields,
        geometry_type,
        crs,
        driver_name,
    )

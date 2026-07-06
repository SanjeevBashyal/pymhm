# -*- coding: utf-8 -*-
"""Standalone compatibility layer for opening the plugin dialog without QGIS.

The real plugin still uses QGIS APIs when QGIS loads it.  The console entry
point installs this module first, so imports such as ``qgis.PyQt`` and
``qgis.gui.QgsMapLayerComboBox`` resolve to small Qt-backed stand-ins.
"""
from __future__ import annotations

import csv
import os
import sys
import types
from enum import IntFlag
from pathlib import Path
from typing import Iterable


_MARKER = "_pymhm_standalone_qgis"
_RASTER_EXTENSIONS = {
    ".asc",
    ".bil",
    ".grd",
    ".img",
    ".nc",
    ".tif",
    ".tiff",
    ".vrt",
}
_VECTOR_EXTENSIONS = {
    ".csv",
    ".dbf",
    ".geojson",
    ".gpkg",
    ".json",
    ".kml",
    ".parquet",
    ".shp",
    ".sqlite",
    ".txt",
}


class _LayerFilter(IntFlag):
    RasterLayer = 1
    VectorLayer = 2
    NoGeometry = 4


def install(force: bool = False) -> bool:
    """Install standalone ``qgis`` and ``processing`` modules when needed.

    Returns ``True`` when the standalone shim is active.  With ``force=False``
    an importable real QGIS installation is left alone.
    """

    existing = sys.modules.get("qgis")
    if existing is not None and getattr(existing, _MARKER, False):
        return True

    if not force:
        try:
            __import__("qgis.core")
            return False
        except Exception:
            pass

    try:
        from PyQt5 import QtCore, QtGui, QtWidgets
    except Exception as exc:  # pragma: no cover - depends on deployment env
        raise RuntimeError(
            "The standalone pymhm GUI requires PyQt5. Install it with "
            "`pip install pymhm[gui]` or install PyQt5 in this environment."
        ) from exc

    if not hasattr(QtCore, "NULL"):
        QtCore.NULL = None

    qgis_mod = types.ModuleType("qgis")
    qgis_mod.__path__ = []
    setattr(qgis_mod, _MARKER, True)

    pyqt_mod = types.ModuleType("qgis.PyQt")
    pyqt_mod.QtCore = QtCore
    pyqt_mod.QtGui = QtGui
    pyqt_mod.QtWidgets = QtWidgets

    core_mod = types.ModuleType("qgis.core")
    gui_mod = types.ModuleType("qgis.gui")

    class QgsRectangle:
        def __init__(self, xmin=0.0, ymin=0.0, xmax=0.0, ymax=0.0):
            self._xmin = float(xmin)
            self._ymin = float(ymin)
            self._xmax = float(xmax)
            self._ymax = float(ymax)

        def xMinimum(self):
            return self._xmin

        def xMaximum(self):
            return self._xmax

        def yMinimum(self):
            return self._ymin

        def yMaximum(self):
            return self._ymax

    class QgsPointXY:
        def __init__(self, x=0.0, y=0.0):
            self._x = float(x)
            self._y = float(y)

        def x(self):
            return self._x

        def y(self):
            return self._y

    class QgsCoordinateReferenceSystem:
        def __init__(self, authid: str | int | None = ""):
            self._authid = ""
            if authid:
                text = str(authid)
                self._authid = text.upper() if text.lower().startswith("epsg:") else text

        def isValid(self):
            return bool(self._authid)

        def authid(self):
            return self._authid

        def postgisSrid(self):
            if self._authid.upper().startswith("EPSG:"):
                try:
                    return int(self._authid.split(":", 1)[1])
                except Exception:
                    return 0
            return 0

        def isGeographic(self):
            return self._authid.upper() in {"EPSG:4326", "EPSG:4258", "EPSG:4269", "CRS:84"}

        def mapUnits(self):
            return "degrees" if self.isGeographic() else "meters"

        def toWkt(self):
            return self._authid

    class QgsCoordinateTransform:
        def __init__(self, source_crs=None, target_crs=None, project=None):
            self.source_crs = source_crs
            self.target_crs = target_crs
            self.project = project

        def setBallparkTransformsAreAppropriate(self, value):
            return None

        def transform(self, point):
            return point

        def transformBoundingBox(self, rectangle):
            return rectangle

    class QgsUnitTypes:
        @staticmethod
        def toAbbreviatedString(unit):
            text = str(unit or "").lower()
            if "degree" in text:
                return "deg"
            if "meter" in text or text == "m":
                return "m"
            return str(unit or "")

    class QgsMapLayer:
        RasterLayer = 0
        VectorLayer = 1

    class QgsMapLayerProxyModel:
        RasterLayer = _LayerFilter.RasterLayer
        VectorLayer = _LayerFilter.VectorLayer
        NoGeometry = _LayerFilter.NoGeometry

    class Qgis:
        LayerFilter = _LayerFilter

    class QgsField:
        def __init__(self, name, field_type=None):
            self._name = str(name)
            self.type = field_type

        def name(self):
            return self._name

    class QgsFields(list):
        def names(self):
            return [field.name() for field in self]

    class QgsGeometry:
        @staticmethod
        def fromPointXY(point):
            geom = QgsGeometry()
            geom._point = point
            return geom

        @staticmethod
        def fromPolylineXY(points):
            geom = QgsGeometry()
            geom._points = list(points)
            return geom

        def asPoint(self):
            return getattr(self, "_point", QgsPointXY())

        def buffer(self, distance, segments):
            return self

        def boundingBox(self):
            point = self.asPoint()
            return QgsRectangle(point.x(), point.y(), point.x(), point.y())

        def distance(self, other):
            return 0.0

        def closestSegmentWithContext(self, point):
            return (0.0, point, 0, 0.0)

    class QgsFeature:
        _next_id = 1

        def __init__(self, fields=None, attributes=None):
            self._id = QgsFeature._next_id
            QgsFeature._next_id += 1
            self._fields = fields or QgsFields()
            self._attributes = {}
            self._geometry = QgsGeometry()
            if isinstance(attributes, dict):
                self._attributes.update(attributes)
            elif attributes is not None:
                self.setAttributes(attributes)

        def id(self):
            return self._id

        def fields(self):
            return self._fields

        def setGeometry(self, geometry):
            self._geometry = geometry

        def geometry(self):
            return self._geometry

        def setAttributes(self, values):
            for index, value in enumerate(values):
                self.setAttribute(index, value)

        def attributes(self):
            names = self._fields.names() if hasattr(self._fields, "names") else []
            return [self._attributes.get(name) for name in names]

        def setAttribute(self, key, value):
            if isinstance(key, int):
                names = self._fields.names() if hasattr(self._fields, "names") else []
                key = names[key] if 0 <= key < len(names) else str(key)
            self._attributes[str(key)] = value

        def attribute(self, key):
            if isinstance(key, int):
                names = self._fields.names() if hasattr(self._fields, "names") else []
                key = names[key] if 0 <= key < len(names) else str(key)
            return self._attributes.get(str(key))

        def __getitem__(self, key):
            return self.attribute(key)

    def _source_path(source: str) -> str:
        text = str(source or "").strip()
        if text.startswith("NETCDF:"):
            first_quote = text.find('"')
            second_quote = text.find('"', first_quote + 1)
            if first_quote >= 0 and second_quote > first_quote:
                return text[first_quote + 1:second_quote]
        return text.split("|", 1)[0].split("?", 1)[0]

    def _csv_rows(path: str):
        try:
            with open(path, "r", encoding="utf-8-sig", newline="") as handle:
                sample = handle.read(4096)
                handle.seek(0)
                try:
                    dialect = csv.Sniffer().sniff(sample)
                except Exception:
                    dialect = csv.excel
                yield from csv.DictReader(handle, dialect=dialect)
        except Exception:
            return

    class _StandaloneLayer:
        _layer_type = QgsMapLayer.VectorLayer

        def __init__(self, source="", name="", provider=""):
            self._source = str(source or "")
            self._name = str(name or "") or os.path.basename(_source_path(self._source))
            self._provider = provider
            self._crs = QgsCoordinateReferenceSystem("EPSG:4326")
            self._fields = None

        def isValid(self):
            path = _source_path(self._source)
            return bool(path) and os.path.exists(path)

        def source(self):
            return self._source

        def name(self):
            return self._name

        def type(self):
            return self._layer_type

        def crs(self):
            return self._crs

        def setCrs(self, crs):
            self._crs = crs

        def extent(self):
            return QgsRectangle()

    class QgsRasterLayer(_StandaloneLayer):
        _layer_type = QgsMapLayer.RasterLayer

        def width(self):
            header = self._ascii_header()
            return int(header.get("ncols", 0))

        def height(self):
            header = self._ascii_header()
            return int(header.get("nrows", 0))

        def extent(self):
            header = self._ascii_header()
            if not header:
                return QgsRectangle()
            xmin = float(header.get("xllcorner", 0.0))
            ymin = float(header.get("yllcorner", 0.0))
            cellsize = float(header.get("cellsize", 0.0))
            xmax = xmin + int(header.get("ncols", 0)) * cellsize
            ymax = ymin + int(header.get("nrows", 0)) * cellsize
            return QgsRectangle(xmin, ymin, xmax, ymax)

        def _ascii_header(self):
            path = _source_path(self._source)
            if not path.lower().endswith(".asc"):
                return {}
            header = {}
            try:
                with open(path, "r", encoding="utf-8") as handle:
                    for _ in range(6):
                        parts = handle.readline().strip().split()
                        if len(parts) >= 2:
                            key = parts[0].lower()
                            value = parts[1]
                            header[key] = float(value)
                return header
            except Exception:
                return {}

    class QgsVectorLayer(_StandaloneLayer):
        _layer_type = QgsMapLayer.VectorLayer

        def fields(self):
            if self._fields is not None:
                return self._fields
            self._fields = QgsFields()
            path = _source_path(self._source)
            rows = _csv_rows(path)
            first = next(rows, None) if rows is not None else None
            if first:
                for name in first.keys():
                    self._fields.append(QgsField(name))
            return self._fields

        def getFeatures(self, feature_ids=None):
            path = _source_path(self._source)
            names = self.fields().names()
            wanted = set(feature_ids or [])
            for row in _csv_rows(path) or []:
                feature = QgsFeature(self.fields(), row)
                if wanted and feature.id() not in wanted:
                    continue
                for name in names:
                    feature.setAttribute(name, row.get(name))
                yield feature

        def getFeature(self, feature_id):
            for feature in self.getFeatures([feature_id]):
                return feature
            return QgsFeature(self.fields())

        def featureCount(self):
            return sum(1 for _ in self.getFeatures())

    class QgsSpatialIndex:
        def __init__(self, features: Iterable | None = None):
            self._features = list(features or [])

        def intersects(self, rectangle):
            return [feature.id() for feature in self._features]

        def nearestNeighbor(self, point, count):
            return [feature.id() for feature in self._features[:count]]

    class QgsWkbTypes:
        Point = 1
        LineString = 2

    class QgsVectorFileWriter:
        NoError = 0

        class SaveVectorOptions:
            def __init__(self):
                self.driverName = "ESRI Shapefile"
                self.fileEncoding = "UTF-8"

        def __init__(self, path, *args, **kwargs):
            self.path = str(path)
            self._error = self.NoError
            self._message = ""
            try:
                Path(self.path).parent.mkdir(parents=True, exist_ok=True)
                Path(self.path).touch()
            except Exception as exc:
                self._error = 1
                self._message = str(exc)

        @classmethod
        def create(cls, path, fields, geometry_type, crs, transform_context, options):
            return cls(path)

        def hasError(self):
            return self._error != self.NoError

        def errorMessage(self):
            return self._message

        def addFeature(self, feature):
            return True

    class QgsProject:
        _instance = None

        def __init__(self):
            self._crs = QgsCoordinateReferenceSystem("EPSG:4326")
            self._layers = {}

        @classmethod
        def instance(cls):
            if cls._instance is None:
                cls._instance = cls()
            return cls._instance

        def crs(self):
            return self._crs

        def setCrs(self, crs):
            self._crs = crs

        def addMapLayer(self, layer):
            self._layers[str(id(layer))] = layer
            return layer

        def mapLayers(self):
            return dict(self._layers)

        def removeMapLayer(self, layer_id):
            self._layers.pop(str(layer_id), None)

        def transformContext(self):
            return None

    class QgsApplication:
        @staticmethod
        def instance():
            return QtWidgets.QApplication.instance()

        @staticmethod
        def processEvents():
            app = QtWidgets.QApplication.instance()
            if app is not None:
                app.processEvents()

    class QgsProjectionSelectionWidget(QtWidgets.QWidget):
        crsChanged = QtCore.pyqtSignal(object)

        def __init__(self, parent=None):
            super().__init__(parent)
            self._crs = QgsCoordinateReferenceSystem("EPSG:4326")
            layout = QtWidgets.QHBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)
            self._combo = QtWidgets.QComboBox(self)
            self._combo.setEditable(True)
            self._combo.addItems(["EPSG:4326", "EPSG:3035", "EPSG:3857"])
            self._combo.setCurrentText(self._crs.authid())
            self._combo.currentTextChanged.connect(self._text_changed)
            layout.addWidget(self._combo)

        def _text_changed(self, text):
            crs = QgsCoordinateReferenceSystem(text.strip())
            if crs.isValid():
                self._crs = crs
                self.crsChanged.emit(crs)

        def setCrs(self, crs):
            if crs and crs.isValid():
                self._crs = crs
                self._combo.setCurrentText(crs.authid())
                self.crsChanged.emit(crs)

        def crs(self):
            return self._crs

    class QgsMapLayerComboBox(QtWidgets.QWidget):
        layerChanged = QtCore.pyqtSignal(object)
        currentIndexChanged = QtCore.pyqtSignal(int)
        currentTextChanged = QtCore.pyqtSignal(str)

        def __init__(self, parent=None):
            super().__init__(parent)
            self._filters = _LayerFilter.RasterLayer | _LayerFilter.VectorLayer
            self._allow_empty = False
            layout = QtWidgets.QHBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(4)
            self._combo = QtWidgets.QComboBox(self)
            self._combo.setEditable(True)
            self._combo.currentIndexChanged.connect(self._index_changed)
            self._combo.currentTextChanged.connect(self._text_changed)
            if self._combo.lineEdit() is not None:
                self._combo.lineEdit().editingFinished.connect(self._emit_layer_changed)
            self._browse = QtWidgets.QToolButton(self)
            self._browse.setText("...")
            self._browse.setToolTip("Browse")
            self._browse.clicked.connect(self._browse_path)
            layout.addWidget(self._combo, 1)
            layout.addWidget(self._browse)

        def setFilters(self, filters):
            try:
                self._filters = _LayerFilter(filters)
            except Exception:
                self._filters = filters

        def setAllowEmptyLayer(self, allow_empty, text=""):
            self._allow_empty = bool(allow_empty)
            if allow_empty and text and self._combo.findText(text) < 0:
                self._combo.insertItem(0, text)

        def setLayer(self, layer):
            if layer is None:
                self.setCurrentIndex(-1)
                self._emit_layer_changed()
                return
            self.setCurrentText(layer.source())
            self._emit_layer_changed(layer)

        def currentLayer(self):
            text = self.currentText().strip()
            if not text:
                return None
            layer_class = self._layer_class(text)
            layer = layer_class(text, os.path.basename(_source_path(text)))
            if not layer.isValid():
                return None
            QgsProject.instance().addMapLayer(layer)
            return layer

        def clear(self):
            self._combo.clear()

        def addItem(self, text, userData=None):
            self._combo.addItem(text, userData)

        def addItems(self, texts):
            self._combo.addItems(texts)

        def count(self):
            return self._combo.count()

        def currentText(self):
            return self._combo.currentText()

        def setCurrentText(self, text):
            text = str(text or "")
            if text and self._combo.findText(text) < 0:
                self._combo.addItem(text)
            self._combo.setCurrentText(text)

        def currentIndex(self):
            return self._combo.currentIndex()

        def setCurrentIndex(self, index):
            self._combo.setCurrentIndex(index)

        def currentData(self):
            return self._combo.currentData()

        def setPlaceholderText(self, text):
            line_edit = self._combo.lineEdit()
            if line_edit is not None:
                line_edit.setPlaceholderText(text)

        def setEditable(self, editable):
            self._combo.setEditable(editable)

        def _index_changed(self, index):
            self.currentIndexChanged.emit(index)
            self._emit_layer_changed()

        def _text_changed(self, text):
            self.currentTextChanged.emit(text)

        def _browse_path(self):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self,
                "Select input file",
                "",
                self._file_filter(),
            )
            if path:
                self.setCurrentText(path)
                self._emit_layer_changed()

        def _emit_layer_changed(self, layer=None):
            if layer is None:
                layer = self.currentLayer()
            self.layerChanged.emit(layer)

        def _layer_class(self, source):
            suffix = Path(_source_path(source)).suffix.lower()
            if suffix in _VECTOR_EXTENSIONS:
                return QgsVectorLayer
            if suffix in _RASTER_EXTENSIONS:
                return QgsRasterLayer
            try:
                if self._filters == _LayerFilter.VectorLayer:
                    return QgsVectorLayer
                if self._filters == _LayerFilter.RasterLayer:
                    return QgsRasterLayer
            except Exception:
                pass
            return QgsVectorLayer

        def _file_filter(self):
            try:
                raster = bool(self._filters & _LayerFilter.RasterLayer)
                vector = bool(self._filters & _LayerFilter.VectorLayer)
            except Exception:
                raster = vector = True
            if raster and not vector:
                return "Raster files (*.tif *.tiff *.asc *.nc *.vrt *.img);;All files (*)"
            if vector and not raster:
                return "Vector/table files (*.shp *.gpkg *.geojson *.csv *.txt *.dbf);;All files (*)"
            return "Supported files (*.tif *.tiff *.asc *.nc *.shp *.gpkg *.geojson *.csv *.txt *.dbf);;All files (*)"

    def _processing_run(name, params):
        raise RuntimeError(
            "QGIS Processing is not available in pymhm standalone mode yet. "
            f"The requested algorithm was '{name}'."
        )

    processing_mod = types.ModuleType("processing")
    processing_mod.run = _processing_run

    for name, value in {
        "Qgis": Qgis,
        "QgsApplication": QgsApplication,
        "QgsCoordinateReferenceSystem": QgsCoordinateReferenceSystem,
        "QgsCoordinateTransform": QgsCoordinateTransform,
        "QgsFeature": QgsFeature,
        "QgsField": QgsField,
        "QgsFields": QgsFields,
        "QgsGeometry": QgsGeometry,
        "QgsMapLayer": QgsMapLayer,
        "QgsMapLayerProxyModel": QgsMapLayerProxyModel,
        "QgsPointXY": QgsPointXY,
        "QgsProject": QgsProject,
        "QgsRasterLayer": QgsRasterLayer,
        "QgsRectangle": QgsRectangle,
        "QgsSpatialIndex": QgsSpatialIndex,
        "QgsUnitTypes": QgsUnitTypes,
        "QgsVectorFileWriter": QgsVectorFileWriter,
        "QgsVectorLayer": QgsVectorLayer,
        "QgsWkbTypes": QgsWkbTypes,
    }.items():
        setattr(core_mod, name, value)

    gui_mod.QgsMapLayerComboBox = QgsMapLayerComboBox
    gui_mod.QgsProjectionSelectionWidget = QgsProjectionSelectionWidget

    qgis_mod.PyQt = pyqt_mod
    qgis_mod.core = core_mod
    qgis_mod.gui = gui_mod

    sys.modules["qgis"] = qgis_mod
    sys.modules["qgis.PyQt"] = pyqt_mod
    sys.modules["qgis.PyQt.QtCore"] = QtCore
    sys.modules["qgis.PyQt.QtGui"] = QtGui
    sys.modules["qgis.PyQt.QtWidgets"] = QtWidgets
    sys.modules["qgis.core"] = core_mod
    sys.modules["qgis.gui"] = gui_mod
    sys.modules["qgsmaplayercombobox"] = gui_mod
    sys.modules["processing"] = processing_mod

    return True


def is_active() -> bool:
    """Return whether the standalone QGIS shim is installed."""

    return bool(getattr(sys.modules.get("qgis"), _MARKER, False))

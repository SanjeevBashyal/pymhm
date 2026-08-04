"""Input-combo helpers for QGIS layers and project files."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

try:
    from qgis.core import (QgsMapLayer, QgsProject, QgsRasterLayer,
                           QgsVectorLayer)
    from qgis.PyQt import QtCore, QtWidgets
except ImportError:
    from .standalone_qgis import install

    install(force=True)
    from qgis.core import (QgsMapLayer, QgsProject, QgsRasterLayer,
                           QgsVectorLayer)
    from qgis.PyQt import QtCore, QtWidgets

from .project_layout import WORKSPACE_FOLDER_NAME
from .pyui.ui_lookup_config_dialog import Ui_LookupConfigDialog

EXCLUDED_PROJECT_FOLDERS = {WORKSPACE_FOLDER_NAME}
INPUT_EXTENSIONS = {
    "dem": {".asc", ".nc", ".tif"},
    "pour_points": {".shp"},
    "land_cover": {".asc", ".nc", ".shp", ".tif"},
    "soil": {".asc", ".nc", ".shp", ".tif"},
    "geology": {".asc", ".nc", ".shp", ".tif"},
    "lai": {".asc", ".nc", ".shp", ".tif"},
    "lookup": {".csv", ".txt"},
}
_KIND_ALIASES = {
    "pour_point": "pour_points",
    "landcover": "land_cover",
    "lc": "land_cover",
}


@dataclass(frozen=True)
class InputItem:
    """One display label and its QComboBox item metadata."""

    label: str
    data: dict


@dataclass(frozen=True)
class LookupConfig:
    """Selected lookup table and its two required fields."""

    lookup_table: str
    mapping_field: str
    class_field: str


def _kind(value: str) -> str:
    kind = _KIND_ALIASES.get(str(value).strip().lower(), str(value).strip().lower())
    if kind not in INPUT_EXTENSIONS:
        raise ValueError(f"Unknown input kind: {value!r}")
    return kind


def scan_project_inputs(project_folder, kind: str) -> list[InputItem]:
    """List eligible project files while pruning plugin-owned top-level folders."""
    project_root = Path(project_folder).expanduser().resolve()
    if not project_root.is_dir():
        return []

    kind = _kind(kind)
    extensions = INPUT_EXTENSIONS[kind]
    items = []
    for current_root, directories, filenames in os.walk(project_root):
        if Path(current_root) == project_root:
            directories[:] = [
                name
                for name in directories
                if name not in EXCLUDED_PROJECT_FOLDERS
            ]
        for filename in filenames:
            path = Path(current_root, filename)
            if path.suffix.lower() not in extensions:
                continue
            label = path.relative_to(project_root).as_posix()
            items.append(
                InputItem(
                    label,
                    {
                        "origin": "file",
                        "kind": kind,
                        "path": str(path.resolve()),
                    },
                )
            )
    return sorted(items, key=lambda item: item.label.casefold())


def scan_project_folders(project_folder) -> list[InputItem]:
    """List project subfolders while pruning the plugin-owned workspace."""
    project_root = Path(project_folder).expanduser().resolve()
    if not project_root.is_dir():
        return []

    items = []
    for current_root, directories, _ in os.walk(project_root):
        current_path = Path(current_root)
        if current_path == project_root:
            excluded = {name.casefold() for name in EXCLUDED_PROJECT_FOLDERS}
            directories[:] = [
                name for name in directories
                if name.casefold() not in excluded
            ]
        for dirname in directories:
            path = current_path / dirname
            items.append(
                InputItem(
                    path.relative_to(project_root).as_posix(),
                    {
                        "origin": "folder",
                        "path": str(path.resolve()),
                    },
                )
            )
    return sorted(items, key=lambda item: item.label.casefold())


def _layer_id(layer):
    getter = getattr(layer, "id", None)
    if not callable(getter):
        return None
    try:
        return getter()
    except Exception:
        return None


def _layer_type(layer):
    getter = getattr(layer, "type", None)
    try:
        return getter() if callable(getter) else None
    except Exception:
        return None


def _accepts_layer(kind: str, layer) -> bool:
    layer_type = _layer_type(layer)
    if kind == "dem":
        return layer_type == QgsMapLayer.RasterLayer
    if kind in {"pour_points", "lookup"}:
        return layer_type == QgsMapLayer.VectorLayer
    return layer_type in {QgsMapLayer.RasterLayer, QgsMapLayer.VectorLayer}


def qgis_layer_item(layer, kind: str) -> InputItem | None:
    """Build an item labelled ``name [CRS]`` for an eligible QGIS layer."""
    kind = _kind(kind)
    if layer is None or not _accepts_layer(kind, layer):
        return None
    try:
        if hasattr(layer, "isValid") and not layer.isValid():
            return None
    except Exception:
        return None

    try:
        name = layer.name() or "Unnamed layer"
    except Exception:
        name = "Unnamed layer"
    authid = ""
    try:
        crs = layer.crs()
        if crs is not None and crs.isValid():
            authid = crs.authid()
    except Exception:
        pass
    label = f"{name} [{authid or 'unknown CRS'}]"
    try:
        source = layer.source()
    except Exception:
        source = ""
    return InputItem(
        label,
        {
            "origin": "qgis",
            "kind": kind,
            "layer": layer,
            "layer_id": _layer_id(layer),
            "source": source,
        },
    )


def loaded_qgis_items(kind: str, project=None) -> list[InputItem]:
    """Return eligible layers currently loaded in a QGIS project."""
    project = project or QgsProject.instance()
    try:
        layers: Iterable = project.mapLayers().values()
    except Exception:
        layers = ()
    items = [qgis_layer_item(layer, kind) for layer in layers]
    return sorted(
        (item for item in items if item is not None),
        key=lambda item: item.label.casefold(),
    )


def _project_layer(layer_id):
    if not layer_id:
        return None
    project = QgsProject.instance()
    getter = getattr(project, "mapLayer", None)
    if callable(getter):
        try:
            layer = getter(layer_id)
            if layer is not None:
                return layer
        except Exception:
            pass
    try:
        for layer in project.mapLayers().values():
            if _layer_id(layer) == layer_id:
                return layer
    except Exception:
        pass
    return None


def _file_layer(path: str, kind: str):
    path = str(path or "")
    if not path:
        return None
    name = Path(path).name
    if kind in {"pour_points", "lookup"} or Path(path).suffix.lower() == ".shp":
        layer = QgsVectorLayer(path, name, "ogr")
    else:
        layer = QgsRasterLayer(path, name)
    try:
        return layer if layer.isValid() else None
    except Exception:
        return None


class InputComboAdapter(QtCore.QObject):
    """Present a plain QComboBox through the legacy layer-combo interface."""

    layerChanged = QtCore.pyqtSignal(object)

    def __init__(self, combo_box, kind: str, parent=None):
        super().__init__(parent or combo_box)
        self.combo_box = combo_box
        self.kind = _kind(kind)
        self.filters = None
        self.allow_empty = True
        self.combo_box.currentIndexChanged.connect(self._emit_layer_changed)

    @property
    def currentIndexChanged(self):
        return self.combo_box.currentIndexChanged

    @property
    def currentTextChanged(self):
        return self.combo_box.currentTextChanged

    def _emit_layer_changed(self, _index=None):
        self.layerChanged.emit(self.currentLayer())

    def currentLayer(self):
        """Return the selected loaded layer or an unregistered file-backed layer."""
        data = self.combo_box.currentData()
        if isinstance(data, dict):
            if data.get("origin") == "qgis":
                layer = _project_layer(data.get("layer_id")) or data.get("layer")
                return layer if layer is not None else None
            if data.get("origin") == "file":
                return _file_layer(data.get("path"), self.kind)
        if data is not None and hasattr(data, "type"):
            return data
        text = self.combo_box.currentText().strip()
        if os.path.isfile(text):
            return _file_layer(text, self.kind)
        return None

    def source_path(self):
        """Return the selected file path or loaded layer provider source."""
        data = self.combo_box.currentData()
        if isinstance(data, dict):
            if data.get("origin") == "file":
                path = data.get("path")
                return str(Path(path).expanduser().resolve()) if path else ""
            if data.get("origin") == "qgis":
                layer = self.currentLayer()
                source = getattr(layer, "source", None)
                return str(source() if callable(source) else source or "")
        layer = self.currentLayer()
        source = getattr(layer, "source", None)
        return str(source() if callable(source) else source or "")

    def setLayer(self, layer):
        """Select an existing layer item, adding it when necessary."""
        if layer is None:
            self.combo_box.setCurrentIndex(-1)
            return

        target_id = _layer_id(layer)
        try:
            target_source = layer.source()
        except Exception:
            target_source = ""
        for index in range(self.combo_box.count()):
            data = self.combo_box.itemData(index)
            if not isinstance(data, dict):
                continue
            if data.get("layer") is layer:
                self.combo_box.setCurrentIndex(index)
                return
            if target_id and data.get("layer_id") == target_id:
                self.combo_box.setCurrentIndex(index)
                return
            if target_source and data.get("source") == target_source:
                self.combo_box.setCurrentIndex(index)
                return

        item = qgis_layer_item(layer, self.kind)
        if item is not None:
            self.combo_box.addItem(item.label, item.data)
            self.combo_box.setCurrentIndex(self.combo_box.count() - 1)

    def setFilters(self, filters):
        self.filters = filters

    def setAllowEmptyLayer(self, allow_empty, text=""):
        self.allow_empty = bool(allow_empty)
        if self.allow_empty and text and self.combo_box.findText(text) < 0:
            self.combo_box.insertItem(0, text, {"origin": "empty"})

    def clear(self):
        self.combo_box.clear()

    def addItem(self, text, userData=None):
        self.combo_box.addItem(text, userData)

    def addItems(self, texts):
        self.combo_box.addItems(texts)

    def count(self):
        return self.combo_box.count()

    def currentText(self):
        return self.combo_box.currentText()

    def setCurrentText(self, text):
        self.combo_box.setCurrentText(text)

    def currentIndex(self):
        return self.combo_box.currentIndex()

    def setCurrentIndex(self, index):
        self.combo_box.setCurrentIndex(index)

    def currentData(self):
        return self.combo_box.currentData()

    def findText(self, text):
        return self.combo_box.findText(text)

    def setEnabled(self, enabled):
        self.combo_box.setEnabled(enabled)

    def isEnabled(self):
        return self.combo_box.isEnabled()

    def setEditable(self, editable):
        self.combo_box.setEditable(editable)

    def setPlaceholderText(self, text):
        line_edit = self.combo_box.lineEdit()
        if line_edit is not None:
            line_edit.setPlaceholderText(text)


def read_lookup_fields(lookup_table) -> list[str]:
    """Return fields using the same table reader as mHM-tools formatting."""
    from .mhm_tools_adapter import read_categorical_lookup_table

    table = read_categorical_lookup_table(lookup_table)
    return [str(column) for column in table.columns if str(column) != "geometry"]


class LookupConfigDialog(QtWidgets.QDialog, Ui_LookupConfigDialog):
    """Select a project lookup table, mapping field, and class field."""

    def __init__(self, project_folder, parent=None, initial=None):
        super().__init__(parent)
        self.setupUi(self)
        self._initial = initial

        for item in scan_project_inputs(project_folder, "lookup"):
            self.lookup_table_combo.addItem(item.label, item.data["path"])

        self.lookup_table_combo.currentIndexChanged.connect(self._refresh_fields)
        self.mapping_field_combo.currentIndexChanged.connect(self._update_ok)
        self.class_field_combo.currentIndexChanged.connect(self._update_ok)
        initial_path = self._initial_value("lookup_table")
        if initial_path:
            for index in range(self.lookup_table_combo.count()):
                if self.lookup_table_combo.itemData(index) == str(initial_path):
                    self.lookup_table_combo.setCurrentIndex(index)
                    break
        self._refresh_fields()

    def _initial_value(self, name):
        if isinstance(self._initial, dict):
            return self._initial.get(name)
        return getattr(self._initial, name, None)

    def _refresh_fields(self, _index=None):
        self.mapping_field_combo.clear()
        self.class_field_combo.clear()
        self.error_label.clear()
        path = self.lookup_table_combo.currentData()
        if not path:
            self._update_ok()
            return
        try:
            fields = read_lookup_fields(path)
        except Exception as error:
            self.error_label.setText(str(error))
            self._update_ok()
            return

        self.mapping_field_combo.addItems(fields)
        self.class_field_combo.addItems(fields)
        for combo, name in (
            (self.mapping_field_combo, "mapping_field"),
            (self.class_field_combo, "class_field"),
        ):
            preferred = self._initial_value(name)
            index = combo.findText(str(preferred)) if preferred else -1
            combo.setCurrentIndex(index)
        self._update_ok()

    def _update_ok(self, _index=None):
        button = self.buttons.button(QtWidgets.QDialogButtonBox.Ok)
        if button is not None:
            button.setEnabled(self.selected_config() is not None)

    def selected_config(self) -> LookupConfig | None:
        path = self.lookup_table_combo.currentData()
        mapping_field = self.mapping_field_combo.currentText().strip()
        class_field = self.class_field_combo.currentText().strip()
        if not path or not mapping_field or not class_field:
            return None
        return LookupConfig(str(path), mapping_field, class_field)

    def accept(self):
        if self.selected_config() is not None:
            super().accept()


__all__ = [
    "EXCLUDED_PROJECT_FOLDERS",
    "INPUT_EXTENSIONS",
    "InputComboAdapter",
    "InputItem",
    "LookupConfig",
    "LookupConfigDialog",
    "loaded_qgis_items",
    "qgis_layer_item",
    "read_lookup_fields",
    "scan_project_inputs",
]

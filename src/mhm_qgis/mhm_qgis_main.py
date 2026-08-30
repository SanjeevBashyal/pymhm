# -*- coding: utf-8 -*-
"""Dialog and UI wiring for the mhm_qgis QGIS plugin."""

import json
import os
import sys
import traceback
from pathlib import Path

# QGIS and PyQt imports
try:
    from qgis.core import (QgsApplication, QgsMapLayer, QgsProject,
                           QgsRasterLayer, QgsVectorLayer)
    from qgis.PyQt import QtCore
    from qgis.PyQt.QtWidgets import QComboBox, QDialog, QFileDialog, QMessageBox
except ImportError:
    from .standalone import install

    install(force=True)
    from qgis.core import (QgsApplication, QgsMapLayer, QgsProject,
                           QgsRasterLayer, QgsVectorLayer)
    from qgis.PyQt import QtCore
    from qgis.PyQt.QtWidgets import QComboBox, QDialog, QFileDialog, QMessageBox

# UI class from the compiled .ui file. The generated module imports
# ``resources_rc`` as a top-level module, so expose the packaged resource module
# under that name when importing through PyPI/package paths.
try:
    from . import resources_rc as _resources_rc  # noqa: F401
    sys.modules.setdefault("resources_rc", _resources_rc)
except Exception:
    pass

from .configuration_processor import ConfigurationProcessor
from .grid_resolution import (build_meteo_l2_grid, ceil_cellsize,
                              display_precision_for_unit, format_resolution,
                              header_bounds, header_for_existing_bounds,
                              is_geographic_unit, l0_header_from_l2,
                              load_meteo_grid_metadata,
                              possible_resolutions, raster_resolution_info,
                              read_header_file)
from .input_selection import (
    INPUT_EXTENSIONS,
    InputComboAdapter,
    LaiNetcdfInputDialog,
    MhmReadyInputDialog,
    SingleLayerInputDialog,
    loaded_qgis_items,
    scan_project_folders,
    scan_project_inputs,
)
from .Meteorology import MeteorologyProcessor
from .Meteorology.forcing import (MeteoFolderSpec, TargetGrid,
                                  inspect_meteo_inputs, resolution_in_crs)
from .Meteorology.inspection_cache import inspect_meteo_folder_cached
from .Meteorology.processor import MeteorologyRun
from .Morphology import MorphologyProcessor
from .Morphology.hydrology.outlets import (
    StationIdError,
    outlet_ids_from_layer,
)
from .project_layout import (data_folder, data_raw_folder,
                             ensure_project_structure, geometry_folder,
                             output_folder, restart_folder,
                             workspace_folder, z_temp_folder)
from .morphology_display import DISPLAY_KEYS
from .standalone import is_active as standalone_is_active
from .qt.bindings.main import bind as bind_main_form
from .qgis_compat import map_layer_filters
from .terminal_dialog import ProjectTerminalDialog
from .task_coordinator import TaskCoordinator
from .thread_display_dialog import ThreadDisplayDialog
from .morphology_task_bridge import MorphologyTaskBridge
from .ui.pyui.ui_mhm_qgis_main import Ui_MhmQgisDialog
# Import utility mixin and processors
from .utils import DialogUtils


class MeteorologyWorkflowWorker(QtCore.QObject):
    """Run only the QGIS-free meteorology phase in a worker thread."""

    log_message = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(str, bool, str)

    def __init__(
            self,
            workflow_key,
            workflow_label,
            meteorology_processor,
            meteorology_run):
        super().__init__()
        self.workflow_key = workflow_key
        self.workflow_label = workflow_label
        self.meteorology_processor = meteorology_processor
        self.meteorology_run = meteorology_run
        self._original_meteo_log_message = None

    @QtCore.pyqtSlot()
    def run(self):
        """Execute meteorology and hand morphology back to the GUI thread."""
        self._original_meteo_log_message = self.meteorology_processor.log_message
        self.meteorology_processor.log_message = self.log_message.emit
        try:
            ok = self.meteorology_processor.process_meteo_forcing(
                self.meteorology_run,
                show_dialog=False,
            )
            self.finished.emit(self.workflow_key, ok, "")
        except Exception as exc:
            self.log_message.emit(
                f"\nERROR: {self.workflow_label} worker failed: {exc}"
            )
            self.log_message.emit(f"Traceback: {traceback.format_exc()}")
            self.finished.emit(self.workflow_key, False, str(exc))
        finally:
            self.meteorology_processor.log_message = (
                self._original_meteo_log_message
            )


class MhmQgisDialog(QDialog, Ui_MhmQgisDialog, DialogUtils):
    def __init__(self, parent=None):
        """Constructor."""

        super(MhmQgisDialog, self).__init__(parent)
        self.setupUi(self)
        self._categorical_lookup_configs = {}
        self._categorical_modes = {}
        self._advanced_inputs = {}
        self._land_cover_ready_source = ""
        self._categorical_ready_configs = {}
        self._lai_input_config = {}
        self._domain_definition_mode = ""
        self._input_adapters = {}
        self.configure_input_adapters()
        self.configure_morphology_display()

        # --- Filter map layer combo boxes to show only relevant layer types ---
        self.input_combo("dem").setFilters(map_layer_filters("RasterLayer"))
        self.input_combo("pour_points").setFilters(map_layer_filters("VectorLayer"))

        # Set filters for new layer combo boxes (both vector and raster allowed)
        self.input_combo("soil").setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.input_combo("land_cover").setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.input_combo("geology").setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.input_combo("lai").setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.configure_input_layer_combo_boxes()

        # --- Instance attributes for managing file paths ---
        self.project_folder = None
        self.geometry_folder = None  # Subfolder for geometry outputs
        self.input_state_filename = "mhm_qgis_input_state.json"
        self._loading_input_state = False
        self._suspend_input_state_saves = False
        self._preserve_missing_layer_state = False
        self._grid_l0_info = None
        self._grid_l2_metadata = None
        self._grid_l2_header = None
        self._preferred_l1_resolution = None
        self._preferred_l11_resolution = None
        self._meteo_inspections = {}
        self._pending_meteo_run = None
        self._terminal_dialog = None
        self._thread_display_dialog = None
        self.task_coordinator = TaskCoordinator(self, max_threads=2)
        self._domain_delineator_dialog = None
        self._morphology_workflow_threads = {}
        self._morphology_workflow_workers = {}
        self._workflow_button_default_styles = {}
        self.spinBox_L2ResolutionMultiplier.setMinimum(1)
        self.spinBox_L2ResolutionMultiplier.setMaximum(1_000_000)
        if self.spinBox_L2ResolutionMultiplier.value() < 1:
            self.spinBox_L2ResolutionMultiplier.setValue(1)
        self.pushButton_executeMeteoMorphSetup.setToolTip(
            "Prepare meteorology, crop and mask morphology, create latlon.nc, "
            "then write morphology ASCII files"
        )
        self._capture_workflow_button_default_styles()

        # --- Initialize processors ---
        self.morphology_processor = MorphologyProcessor(self)
        self.meteorology_processor = MeteorologyProcessor(self)
        self.configuration_processor = ConfigurationProcessor(self)
        self.morphology_tasks = MorphologyTaskBridge(self)

        # --- Connect signals and slots ---
        bind_main_form(self)
        self.refresh_input_sources()
        self.refresh_meteo_display()
        self.refresh_meteo_folder_sources()
        self.refresh_grid_resolution_controls()

    def configure_input_adapters(self):
        """Wrap each plain input combo in the layer-combo interface the code uses."""

        input_widgets = {
            "dem": "comboBox_demInput",
            "pour_points": "comboBox_pourPointInput",
            "land_cover": "comboBox_landUseInput",
            "soil": "comboBox_soilInput",
            "geology": "comboBox_geologyInput",
            "lai": "comboBox_laiInput",
        }
        for kind, combo_name in input_widgets.items():
            combo = getattr(self, combo_name, None)
            if combo is None:
                if kind not in {"land_cover", "soil", "geology", "lai"}:
                    continue
                combo = QComboBox(self)
                combo.setObjectName(combo_name)
                combo.hide()
                setattr(self, combo_name, combo)
            self._input_adapters[kind] = InputComboAdapter(combo, kind, self)

    def input_combo(self, kind):
        """Return the layer-combo adapter for one input kind, or None."""
        return self._input_adapters.get(kind)

    def refresh_meteo_display(self):
        """Re-read the prepared forcing and re-range the display controls."""
        from .qgis_bridge.display import meteo as meteo_display

        meteo_display.refresh(self)

    def configure_morphology_display(self):
        """Attach stable keys to the morphology display choices."""
        combo = self.comboBox_morphVariableToDisplay
        if combo is not None:
            for index, key in enumerate(DISPLAY_KEYS[: combo.count()]):
                combo.setItemData(index, key)
        editor = self.dateTimeEdit_forHistoricalInputs
        if editor is not None:
            editor.setEnabled(False)

    def configure_input_layer_combo_boxes(self):
        """Allow input layer boxes to start empty so layers are chosen deliberately."""

        for adapter in self._input_adapters.values():
            adapter.setAllowEmptyLayer(True)
            adapter.setLayer(None)

    def refresh_input_sources(self):
        """Populate input boxes from QGIS and, when enabled, the project folder."""
        include_files = bool(
            self.project_folder and self.checkBox_enableFolderSearch.isChecked()
        )
        for kind, adapter in self._input_adapters.items():
            combo = adapter.combo_box
            previous = combo.currentData()
            combo.blockSignals(True)
            combo.clear()
            for item in loaded_qgis_items(kind):
                combo.addItem(item.label, item.data)
            if include_files:
                for item in scan_project_inputs(self.project_folder, kind):
                    combo.addItem(item.label, item.data)
            index = self._matching_input_index(combo, previous)
            if (
                index < 0
                and isinstance(previous, dict)
                and previous.get("origin") == "file"
                and previous.get("manual")
                and os.path.isfile(previous.get("path", ""))
            ):
                combo.addItem(
                    previous.get("label") or previous["path"],
                    previous,
                )
                index = combo.count() - 1
            combo.setCurrentIndex(index)
            combo.blockSignals(False)
            if isinstance(previous, dict) and index < 0:
                adapter.layerChanged.emit(adapter.currentLayer())

    def populate_pour_point_outlet_fields(self, layer=None, preferred=None):
        """Populate the outlet ID field selector from the pour-point layer."""
        combo = self.comboBox_pourPointOutletID
        previous = str(preferred or combo.currentText() or "").strip()
        if layer is None:
            layer = self.input_combo("pour_points").currentLayer()

        names = []
        try:
            if layer and layer.isValid():
                names = list(layer.fields().names())
        except Exception:
            names = []

        combo.blockSignals(True)
        combo.clear()
        combo.addItem("")
        combo.addItems(names)
        index = combo.findText(previous) if previous else -1
        if index < 0:
            for candidate in range(1, combo.count()):
                if combo.itemText(candidate).casefold() == previous.casefold():
                    index = candidate
                    break
        if index < 0:
            for candidate in range(1, combo.count()):
                if combo.itemText(candidate).casefold() == "station_id":
                    index = candidate
                    break
        combo.setCurrentIndex(index if index >= 0 else 0)
        combo.blockSignals(False)

    def selected_outlet_id_field(self):
        """Return the selected pour-point unique ID field."""
        return self.comboBox_pourPointOutletID.currentText().strip()

    def selected_outlet_ids(self):
        """Return validated unique outlet IDs from the selected field."""
        layer = self.input_combo("pour_points").currentLayer()
        field_name = self.selected_outlet_id_field()
        if not field_name:
            raise StationIdError("Select the pour-point outlet ID field.")
        return outlet_ids_from_layer(layer, field_name)

    def selected_input_file_paths(self):
        """Return local files already selected elsewhere in the plugin."""
        paths = set()
        for _, widget in self.input_layer_widgets():
            source = getattr(widget, "source_path", lambda: "")()
            source = str(source or "").split("|", 1)[0]
            if os.path.isfile(source):
                paths.add(os.path.normcase(os.path.abspath(source)))
        for config in self._categorical_lookup_configs.values():
            path = str(config.get("lookup_table", "") or "")
            if os.path.isfile(path):
                paths.add(os.path.normcase(os.path.abspath(path)))
        for _, widget in self.input_text_widgets():
            path = widget.text().strip()
            if os.path.isfile(path):
                paths.add(os.path.normcase(os.path.abspath(path)))
        for value in self._advanced_inputs.values():
            data = value.as_dict() if hasattr(value, "as_dict") else value
            for path in self._nested_existing_paths(data):
                paths.add(os.path.normcase(os.path.abspath(path)))
        if os.path.isfile(self._land_cover_ready_source):
            paths.add(
                os.path.normcase(os.path.abspath(self._land_cover_ready_source))
            )
        return paths

    @classmethod
    def _nested_existing_paths(cls, value):
        if isinstance(value, dict):
            for item in value.values():
                yield from cls._nested_existing_paths(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                yield from cls._nested_existing_paths(item)
        elif isinstance(value, str) and os.path.isfile(value):
            yield value

    def meteo_input_widgets(self):
        """Return folder/source widgets for the three meteorology inputs."""
        return (
            (
                "precipitation",
                self.comboBox_precipitationFile,
                self.comboBox_precipitationDataSource,
            ),
            (
                "temperature",
                self.comboBox_temperatureFile,
                self.comboBox_temperatureDataSource,
            ),
            ("pet", self.comboBox_petFile, self.comboBox_petDataSource),
        )

    def refresh_meteo_folder_sources(self):
        """Populate meteo folder boxes from the outer project directory."""
        available = (
            scan_project_folders(self.project_folder)
            if self.project_folder else ()
        )
        for _, combo, _ in self.meteo_input_widgets():
            previous = self.selected_folder_path(combo)
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("", None)
            for item in available:
                combo.addItem(item.label, item.data)
            index = self._folder_combo_index(combo, previous)
            if index < 0 and previous and os.path.isdir(previous):
                self._add_folder_combo_item(combo, previous)
                index = combo.count() - 1
            combo.setCurrentIndex(index if index >= 0 else 0)
            combo.blockSignals(False)

    def selected_meteo_folder(self, kind):
        """Return the selected absolute folder path for one meteo input."""
        for input_kind, combo, _ in self.meteo_input_widgets():
            if input_kind == kind:
                return self.selected_folder_path(combo)
        return ""

    def selected_meteo_source(self, kind):
        """Return a stable internal source token for one meteo input."""
        for input_kind, _, combo in self.meteo_input_widgets():
            if input_kind != kind:
                continue
            text = combo.currentText().strip().lower()
            if text.replace("_", " ") == "mhm ready":
                return "mhm_ready"
            if text.replace("-", "").replace("_", "") == "era5land":
                return "era5land"
        return ""

    @staticmethod
    def selected_folder_path(combo):
        data = combo.currentData()
        if isinstance(data, dict):
            path = data.get("path", "")
        else:
            path = data if isinstance(data, str) else ""
        return os.path.abspath(path) if path else ""

    @staticmethod
    def _folder_combo_index(combo, path):
        if not path:
            return -1
        normalized = os.path.normcase(os.path.abspath(path))
        for index in range(combo.count()):
            data = combo.itemData(index)
            candidate = data.get("path", "") if isinstance(data, dict) else data
            if candidate and os.path.normcase(os.path.abspath(candidate)) == normalized:
                return index
        return -1

    def _add_folder_combo_item(self, combo, folder):
        folder = os.path.abspath(folder)
        label = folder
        if self.project_folder:
            try:
                relative = os.path.relpath(folder, self.project_folder)
                if relative != ".." and not relative.startswith(f"..{os.sep}"):
                    label = relative.replace("\\", "/")
            except ValueError:
                pass
        combo.addItem(label, {"origin": "folder", "path": folder, "manual": True})

    def browse_meteo_input_folder(self, kind):
        """Browse for, add, and select one meteorology input folder."""
        folder = QFileDialog.getExistingDirectory(
            self,
            f"Select {kind.title()} Data Folder",
            self.selected_meteo_folder(kind) or self.project_folder or "",
        )
        if not folder:
            return
        for input_kind, combo, _ in self.meteo_input_widgets():
            if input_kind != kind:
                continue
            index = self._folder_combo_index(combo, folder)
            if index < 0:
                self._add_folder_combo_item(combo, folder)
                index = combo.count() - 1
            combo.setCurrentIndex(index)
            self.log_message(
                f"{kind.title()} data folder selected: {os.path.abspath(folder)}"
            )
            return

    def handle_meteo_input_changed(self, kind):
        """Persist and inspect a changed meteorology folder/source pair."""
        if self._loading_input_state:
            return
        self.save_input_state()
        self.invalidate_meteo_morph_setup()
        self.inspect_meteo_selection(kind, show_errors=True)

    def selected_meteo_crs(self):
        """Return the selected CRS in a form accepted by pyproj."""
        crs = self.get_crs()
        if crs is None or not crs.isValid():
            return ""
        authid = crs.authid()
        if authid:
            return authid
        to_wkt = getattr(crs, "toWkt", None)
        return to_wkt() if callable(to_wkt) else ""

    def meteo_folder_spec(self, kind, required=True):
        """Return one validated folder/source selection."""
        folder = self.selected_meteo_folder(kind)
        source = self.selected_meteo_source(kind)
        if not folder and not required:
            return None
        if not folder:
            raise ValueError(f"Select the {kind} data folder.")
        if not source:
            raise ValueError(f"Select the {kind} data source.")
        if not os.path.isdir(folder):
            raise ValueError(f"The {kind} data folder does not exist:\n{folder}")
        return MeteoFolderSpec(
            kind=kind,
            folder=Path(folder),
            source=source,
            crs=(
                self.selected_meteo_crs()
                if source == "mhm_ready" else None
            ),
        )

    def selected_meteo_specs(self):
        """Return required precipitation/temperature and optional PET inputs."""
        precipitation = self.meteo_folder_spec("precipitation")
        temperature = self.meteo_folder_spec("temperature")
        pet = self.meteo_folder_spec("pet", required=False)
        return precipitation, temperature, pet

    def clear_precipitation_resolution_labels(self):
        """Clear the raw precipitation resolution display."""
        for name in (
                "label_precipitationResolutionValue",
                "label_precipitationResolutionUnit",
                "label_precipitationResolutionMultiplier"):
            label = getattr(self, name, None)
            if label is not None:
                label.setText("")

    def inspect_meteo_selection(self, kind, show_errors=False):
        """Validate one selected meteo folder and refresh its metadata."""
        folder = self.selected_meteo_folder(kind)
        source = self.selected_meteo_source(kind)
        if not folder or not source:
            self._meteo_inspections.pop(kind, None)
            if kind == "precipitation":
                self.clear_precipitation_resolution_labels()
            return None

        try:
            metadata = inspect_meteo_folder_cached(
                self.project_folder,
                folder,
                kind,
                source,
                crs_fallback=(
                    self.selected_meteo_crs()
                    if source == "mhm_ready" else None
                ),
                log=self.log_message,
            )
            self._meteo_inspections[kind] = metadata
            if kind == "precipitation":
                converted = resolution_in_crs(
                    metadata,
                    self.selected_meteo_crs() or None,
                )
                self.label_precipitationResolutionValue.setText(
                    format_resolution(converted.resolution, converted.unit)
                )
                self.label_precipitationResolutionUnit.setText(converted.unit)
                l0_resolution = self.current_l0_resolution()
                multiplier = (
                    f"{converted.resolution / l0_resolution:.1f}"
                    if l0_resolution else ""
                )
                self.label_precipitationResolutionMultiplier.setText(multiplier)
            self.log_message(
                f"{kind.title()} metadata: {len(metadata.files)} NetCDF file(s), "
                f"{metadata.shape[1]} x {metadata.shape[0]} cells."
            )
            return metadata
        except Exception as error:
            self._meteo_inspections.pop(kind, None)
            if kind == "precipitation":
                self.clear_precipitation_resolution_labels()
            self.log_message(f"ERROR: Invalid {kind} meteorology input: {error}")
            if show_errors:
                QMessageBox.warning(
                    self,
                    f"Invalid {kind.title()} Data",
                    str(error),
                )
            return None

    @staticmethod
    def _matching_input_index(combo, previous):
        if not isinstance(previous, dict):
            return -1
        origin = previous.get("origin")
        identity = (
            previous.get("path")
            if origin == "file"
            else previous.get("layer_id") or previous.get("source")
        )
        for index in range(combo.count()):
            data = combo.itemData(index)
            if not isinstance(data, dict) or data.get("origin") != origin:
                continue
            candidate = (
                data.get("path")
                if origin == "file"
                else data.get("layer_id") or data.get("source")
            )
            if candidate == identity:
                return index
        return -1

    def browse_input_file(self, kind):
        """Add a manually selected file to one input box."""
        patterns = " ".join(
            f"*{suffix}" for suffix in sorted(INPUT_EXTENSIONS[kind])
        )
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Select {kind.replace('_', ' ').title()} Input",
            self.project_folder or "",
            f"Supported files ({patterns});;All files (*)",
        )
        if not path:
            return
        path = os.path.abspath(path)
        if os.path.splitext(path)[1].lower() not in INPUT_EXTENSIONS[kind]:
            QMessageBox.warning(
                self,
                "Unsupported Input",
                f"Select a {kind.replace('_', ' ')} file with one of these "
                f"extensions: {', '.join(sorted(INPUT_EXTENSIONS[kind]))}.",
            )
            return
        label = path
        if self.project_folder:
            try:
                relative = os.path.relpath(path, self.project_folder)
                if relative != ".." and not relative.startswith(f"..{os.sep}"):
                    label = relative.replace("\\", "/")
            except ValueError:
                pass
        adapter = self._input_adapters[kind]
        combo = adapter.combo_box
        for index in range(combo.count()):
            data = combo.itemData(index)
            if isinstance(data, dict) and data.get("path") == path:
                combo.setCurrentIndex(index)
                return
        combo.addItem(
            label,
            {
                "origin": "file",
                "kind": kind,
                "path": path,
                "label": label,
                "manual": True,
            },
        )
        combo.setCurrentIndex(combo.count() - 1)

    def connect_optional_processor_button(self, name, label, callback):
        """Connect a processing control when it exists in the active UI."""
        button = getattr(self, name, None)
        if button is not None:
            self.connect_processor_button(
                button, label, callback, background=True
            )

    def handle_l2_multiplier_changed(self, value=None):
        """Invalidate prepared setup state when the requested L2 grid changes."""
        if self._loading_input_state:
            return
        self.save_input_state()
        self.invalidate_meteo_morph_setup()

    def handle_model_input_changed(self, value=None):
        """Persist model inputs and invalidate a previously completed setup."""
        if self._loading_input_state:
            return
        self.save_input_state()
        self.invalidate_meteo_morph_setup()

    def refresh_grid_resolution_controls(self):
        """Refresh L0, L2, L1, and L11 controls from current project state."""
        self.update_l0_resolution_from_dem()
        self.update_l2_resolution_from_metadata()
        self.refresh_l1_l11_resolution_options()

    def update_l0_resolution_from_dem(self, layer=None):
        """Prepare/read the filled DEM and show its resolution as L0."""
        info = self.filled_dem_resolution_info()
        self._grid_l0_info = info
        if not info:
            self._set_resolution_labels("L0", "", "")
            self.refresh_l1_l11_resolution_options()
            self.update_latlon_button_state()
            return

        self._set_resolution_labels(
            "L0",
            format_resolution(info["resolution"], info["unit"]),
            info["unit"],
        )
        if abs(info["x_resolution"] - info["y_resolution"]) > max(info["resolution"], 1.0) * 1e-6:
            self.log_message(
                "WARNING: DEM pixels are not square. L0 uses the average of "
                f"x={info['x_resolution']} and y={info['y_resolution']}.")
        self.refresh_l1_l11_resolution_options()

    def filled_dem_resolution_info(self):
        """Return L0 resolution from the prepared filled DEM raster."""
        if not self.project_folder:
            return None
        if not self.input_combo("dem").currentLayer():
            return None

        processor = getattr(self, "morphology_processor", None)
        if processor is None:
            return None

        try:
            if not processor._ensure_filled_dem(processor.fill_dem):
                return None
        except Exception as e:
            self.log_message(f"WARNING: Could not prepare filled DEM for L0 resolution: {e}")
            return None

        filled_path = getattr(processor, "filled_dem_path", None)
        if not filled_path or not os.path.exists(filled_path):
            filled_path = os.path.join(
                geometry_folder(self.project_folder),
                "1_dem_filled.tif",
            )
            if os.path.exists(filled_path):
                processor.filled_dem_path = filled_path

        if not filled_path or not os.path.exists(filled_path):
            return None

        filled_layer = QgsRasterLayer(filled_path, "Filled_DEM")
        if not filled_layer.isValid():
            self.log_message("WARNING: Filled DEM exists but could not be read for L0 resolution.")
            return None
        return raster_resolution_info(filled_layer)

    def update_l2_resolution_from_metadata(self, metadata=None):
        """Show L2 resolution from saved meteo grid metadata or existing headers."""
        if metadata is None and self.project_folder:
            metadata = load_meteo_grid_metadata(self.project_folder)

        header = None
        unit = ""
        if metadata:
            header = metadata.get("l2_header")
            unit = metadata.get("l2_unit", "")
            if not unit and self._grid_l0_info:
                unit = self._grid_l0_info.get("unit", "")
            if header:
                # Keep the exact L2 cell size. It is n x the filled DEM cell
                # size by construction, and re-rounding it here breaks that
                # relationship for repeating values such as 120 x 1/1200 deg.
                header = dict(header)
                header["unit"] = unit
                metadata["l2_header"] = header
                metadata["l2_resolution"] = header["cellsize"]
                metadata["l2_unit"] = unit
        elif self.project_folder:
            unit = self._grid_l0_info.get("unit", "") if self._grid_l0_info else ""
            header = read_header_file(
                os.path.join(data_folder(self.project_folder), "meteo", "pre", "header.txt"),
                unit=unit,
            )
            if header:
                metadata = {
                    "l2_resolution": header["cellsize"],
                    "l2_unit": unit,
                    "l2_header": header,
                }

        self._grid_l2_metadata = metadata
        self._grid_l2_header = header
        if not metadata:
            self.update_extent_labels()
            self.refresh_l1_l11_resolution_options()
            self.update_latlon_button_state()
            return

        self.update_extent_labels(metadata)
        self.refresh_l1_l11_resolution_options()
        self.update_latlon_button_state()

    def set_meteo_l2_grid_metadata(self, metadata):
        """Store freshly prepared L2 metadata and update resolution controls."""
        self.update_l2_resolution_from_metadata(metadata)
        self.save_input_state()

    def prepare_meteo_l2_grid(self, precipitation_metadata=None):
        """Build the requested L2 grid used by meteorology processing."""
        self.update_l0_resolution_from_dem()
        raw = None
        if precipitation_metadata is not None:
            converted = resolution_in_crs(
                precipitation_metadata,
                self.selected_meteo_crs() or None,
            )
            raw = {
                "resolution": converted.resolution,
                "x_resolution": converted.x_resolution,
                "y_resolution": converted.y_resolution,
                "unit": converted.unit,
                "source_file": str(precipitation_metadata.files[0]),
            }
        grid = build_meteo_l2_grid(
            self,
            self.spinBox_L2ResolutionMultiplier.value(),
            raw_metadata=raw,
        )
        metadata = grid.get("metadata", {})
        self.log_message(
            "Meteo L2 grid: "
            f"{format_resolution(metadata.get('l2_resolution'), metadata.get('l2_unit', ''))} "
            f"{metadata.get('l2_unit', '')} "
            f"({metadata.get('l2_ratio_to_l0')} x L0)."
        )
        return grid

    def refresh_l1_l11_resolution_options(self):
        """Populate L1 and L11 resolution choices from L0/L2 compatibility."""
        l0_resolution = self.current_l0_resolution()
        l2_resolution = self.current_l2_resolution()
        unit = self.current_grid_unit()

        if not self._grid_l2_header or not l2_resolution:
            self.disable_l1_l11_resolution_options()
            self.update_latlon_button_state()
            return

        l1_values = []
        if l0_resolution:
            l1_values = possible_resolutions(l0_resolution, l2_resolution, unit)

        preferred_l1 = self._preferred_l1_resolution
        matched_l1 = self._populate_resolution_combo(
            self.comboBox_L1,
            l1_values,
            preferred_l1,
            unit,
        )
        if preferred_l1 is not None and l1_values and l0_resolution and l2_resolution:
            if not matched_l1:
                self.log_resolution_preference_warning(
                    "L1",
                    preferred_l1,
                    l1_values[0],
                    unit,
                )
            self._preferred_l1_resolution = None
        self.label_L1ResolutionUnit.setText(unit if l1_values else "")

        self.handle_l1_resolution_changed()

    def handle_l1_resolution_changed(self):
        """Refresh L1 label and rebuild L11 choices for the selected L1."""
        self.update_l1_resolution_label()
        l1_resolution = self.current_l1_resolution()
        l2_resolution = self.current_l2_resolution()
        unit = self.current_grid_unit()

        if not self._grid_l2_header or not l2_resolution or not l1_resolution:
            self.disable_l11_resolution_options()
            self.update_latlon_button_state()
            return

        l11_values = []
        l11_values = possible_resolutions(l1_resolution, l2_resolution, unit)

        preferred_l11 = self._preferred_l11_resolution
        matched_l11 = self._populate_resolution_combo(
            self.comboBox_L11,
            l11_values,
            preferred_l11,
            unit,
        )
        if preferred_l11 is not None and l11_values and l1_resolution and l2_resolution:
            if not matched_l11:
                self.log_resolution_preference_warning(
                    "L11",
                    preferred_l11,
                    l11_values[0],
                    unit,
                )
            self._preferred_l11_resolution = None
        self.label_L11ResolutionUnit.setText(unit if l11_values else "")
        self.update_l11_resolution_label()
        self.update_latlon_button_state()

    def disable_l1_l11_resolution_options(self):
        """Disable L1/L11 controls until meteo L2 grid metadata exists."""
        self.disable_resolution_combo(self.comboBox_L1)
        self.disable_l11_resolution_options()
        self._set_resolution_labels("L1", "", "")

    def disable_l11_resolution_options(self):
        """Disable L11 controls and clear its labels."""
        self.disable_resolution_combo(self.comboBox_L11)
        self._set_resolution_labels("L11", "", "")

    def disable_resolution_combo(self, combo_box):
        """Clear and disable a grid-resolution combo box."""
        if combo_box is None:
            return
        try:
            combo_box.blockSignals(True)
        except Exception:
            pass
        combo_box.clear()
        combo_box.setEnabled(False)
        try:
            combo_box.blockSignals(False)
        except Exception:
            pass

    def update_l1_resolution_label(self):
        """Show the selected L1 relation to L0."""
        value = self.current_l1_resolution()
        l0_resolution = self.current_l0_resolution()
        label = self.label_L1Resolution
        if label is None:
            return
        if value and l0_resolution:
            label.setText(f"{value / l0_resolution:g} x L0")
        else:
            label.setText("")

    def update_l11_resolution_label(self):
        """Show the selected L11 relation to L1."""
        value = self.current_l11_resolution()
        l1_resolution = self.current_l1_resolution()
        label = self.label_L11Resolution
        if label is None:
            return
        if value and l1_resolution:
            label.setText(f"{value / l1_resolution:g} x L1")
        else:
            label.setText("")
        self.update_latlon_button_state()

    def current_l0_resolution(self):
        """Return the current L0 cell size, unrounded when it is known."""
        if self._grid_l0_info:
            return float(
                self._grid_l0_info.get("exact_resolution")
                or self._grid_l0_info["resolution"]
            )
        return None

    def current_l2_resolution(self):
        """Return current L2 resolution."""
        if self._grid_l2_metadata and self._grid_l2_metadata.get("l2_resolution"):
            return float(self._grid_l2_metadata["l2_resolution"])
        if self._grid_l2_header:
            return float(self._grid_l2_header["cellsize"])
        return None

    def current_l1_resolution(self):
        """Return selected L1 resolution."""
        return self._current_combo_resolution(self.comboBox_L1)

    def current_l11_resolution(self):
        """Return selected L11 resolution."""
        return self._current_combo_resolution(self.comboBox_L11)

    def current_grid_unit(self):
        """Return the unit shared by the derived grid-resolution controls."""
        if self._grid_l2_metadata and self._grid_l2_metadata.get("l2_unit"):
            return self._grid_l2_metadata["l2_unit"]
        if self._grid_l0_info:
            return self._grid_l0_info.get("unit", "")
        return ""

    def grid_level_headers(self):
        """Return compatible L0, L1, L11, and L2 headers for latlon.nc."""
        if not self._grid_l2_header:
            self.update_l2_resolution_from_metadata()
        if not self._grid_l2_header:
            raise ValueError("L2 grid is not available. Run and save meteorology data first.")

        l0_resolution = self.current_l0_resolution()
        l1_resolution = self.current_l1_resolution()
        l11_resolution = self.current_l11_resolution()
        if not l0_resolution or not l1_resolution or not l11_resolution:
            raise ValueError("L0, L1, and L11 resolutions must be available.")

        unit = self.current_grid_unit()
        l2_header = dict(self._grid_l2_header)
        l2_header["unit"] = unit
        ratio = int(
            (self._grid_l2_metadata or {}).get("l2_ratio_to_l0")
            or round(float(l2_header["cellsize"]) / float(l0_resolution))
        )
        return {
            "L0": l0_header_from_l2(l2_header, l0_resolution, ratio),
            "L1": header_for_existing_bounds(l2_header, l1_resolution, unit),
            "L11": header_for_existing_bounds(l2_header, l11_resolution, unit),
            "L2": l2_header,
        }

    def grid_configuration_snapshot(self):
        """Return current grid configuration for project-state reference."""
        snapshot = {
            "l0_resolution": self.current_l0_resolution(),
            "l1_resolution": self.current_l1_resolution(),
            "l11_resolution": self.current_l11_resolution(),
            "l2_resolution": self.current_l2_resolution(),
            "unit": self.current_grid_unit(),
        }
        try:
            snapshot["headers"] = self.grid_level_headers()
        except Exception:
            snapshot["headers"] = {}
        return snapshot

    def update_extent_labels(self, metadata=None):
        """Show the final model extent from the prepared L2 grid header."""
        label_names = (
            "label_minimumEasting",
            "label_maximumEasting",
            "label_minimumNorthing",
            "label_maximumNorthing",
        )
        labels = [getattr(self, name, None) for name in label_names]
        if any(label is None for label in labels):
            return

        header = None
        if metadata:
            header = metadata.get("l2_header")
        if header is None:
            header = self._grid_l2_header

        if not header:
            for label in labels:
                label.setText("")
            return

        unit = (
            (metadata or {}).get("l2_unit")
            or (self._grid_l2_metadata or {}).get("l2_unit")
            or self.current_grid_unit()
        )
        precision = display_precision_for_unit(unit)
        xmin, xmax, ymin, ymax = header_bounds(header)
        values = (xmin, xmax, ymin, ymax)
        for label, value in zip(labels, values):
            label.setText(format_resolution(value, unit, precision=precision))

    def update_latlon_button_state(self):
        """Enable latlon creation once all derived grid levels are available."""
        button = self.pushButton_createLatLon
        if button is None:
            return

        enabled = bool(
            self.current_l0_resolution()
            and self.current_l1_resolution()
            and self.current_l11_resolution()
            and self._grid_l2_header
        )
        button.setEnabled(enabled)

    def _set_resolution_labels(self, level, value_text, unit_text):
        """Set value and unit labels for a grid level."""
        value_label = getattr(self, f"label_{level}Resolution", None)
        unit_label = getattr(self, f"label_{level}ResolutionUnit", None)
        if value_label is not None:
            value_label.setText(value_text or "")
        if unit_label is not None:
            unit_label.setText(unit_text or "")

    def _populate_resolution_combo(self, combo_box, values, preferred_value=None, unit=None):
        """Populate a resolution combo box while preserving a compatible selection."""
        if combo_box is None:
            return False
        current_value = self._current_combo_resolution(combo_box)
        preferred = preferred_value or current_value
        try:
            combo_box.blockSignals(True)
        except Exception:
            pass
        combo_box.clear()
        for value in values:
            combo_box.addItem(format_resolution(value, unit), float(value))
        matched_preferred = preferred is None
        if values:
            selected_index = 0
            if preferred:
                for index, value in enumerate(values):
                    if self.resolution_values_match(value, preferred, unit):
                        selected_index = index
                        matched_preferred = True
                        break
            combo_box.setCurrentIndex(selected_index)
        combo_box.setEnabled(True)
        try:
            combo_box.blockSignals(False)
        except Exception:
            pass
        return matched_preferred

    def resolution_state_precision(self, unit):
        """Return precision used when restoring saved grid-resolution choices."""
        if is_geographic_unit(unit):
            return 8
        return 6

    def normalized_state_resolution(self, value, unit):
        """Normalize saved and available resolutions for robust state matching."""
        return ceil_cellsize(value, unit, precision=self.resolution_state_precision(unit))

    def resolution_values_match(self, available_value, preferred_value, unit):
        """Return True when two resolutions are equivalent for state restore."""
        try:
            available = self.normalized_state_resolution(available_value, unit)
            preferred = self.normalized_state_resolution(preferred_value, unit)
        except (TypeError, ValueError):
            return False
        tolerance = 10 ** (-self.resolution_state_precision(unit))
        return abs(available - preferred) <= tolerance

    def log_resolution_preference_warning(
            self,
            level,
            saved_value,
            fallback_value,
            unit):
        """Log when a saved L1/L11 choice is no longer selectable."""
        precision = self.resolution_state_precision(unit)
        self.log_message(
            "WARNING: Saved "
            f"{level} resolution {format_resolution(saved_value, unit, precision)} "
            f"{unit or ''} is not compatible with the current grid. "
            f"Using {format_resolution(fallback_value, unit)} {unit or ''}."
        )

    def _current_combo_resolution(self, combo_box):
        """Return current numeric resolution from a combo box."""
        if combo_box is None or combo_box.count() == 0:
            return None
        try:
            data = combo_box.currentData()
            if data is not None:
                return float(data)
        except Exception:
            pass
        try:
            return float(combo_box.currentText())
        except (TypeError, ValueError):
            return None

    def connect_input_state_teardown_guards(self):
        """Avoid overwriting saved layer selections during QGIS teardown."""
        application = QtCore.QCoreApplication.instance()
        if application is not None:
            try:
                application.aboutToQuit.connect(self.suspend_input_state_saves)
            except Exception:
                pass

        try:
            project = QgsProject.instance()
        except Exception:
            project = None
        if project is None:
            return

        removal_start_signal = getattr(project, "layersWillBeRemoved", None)
        if removal_start_signal is not None:
            try:
                removal_start_signal.connect(self.preserve_missing_layer_state)
            except Exception:
                pass

        for signal_name in ("layersRemoved", "cleared"):
            signal = getattr(project, signal_name, None)
            if signal is None:
                continue
            try:
                signal.connect(self.release_missing_layer_state_preservation)
            except Exception:
                pass

    def suspend_input_state_saves(self, *args):
        """Stop input-state writes once QGIS is shutting down."""
        self._preserve_missing_layer_state = True
        self._suspend_input_state_saves = True

    def preserve_missing_layer_state(self, *args):
        """Keep saved layer entries if QGIS removes layers from combo boxes."""
        self._preserve_missing_layer_state = True

    def release_missing_layer_state_preservation(self, *args):
        """Stop preserving empty layer selections after a removal batch ends."""
        if self._suspend_input_state_saves:
            return
        QtCore.QTimer.singleShot(0, self.clear_missing_layer_state_preservation)

    def clear_missing_layer_state_preservation(self):
        """Return to normal layer-state saving after project layer removal."""
        if not self._suspend_input_state_saves:
            self._preserve_missing_layer_state = False

    def connect_processor_button(
        self, button, action_name, callback, *, background=False
    ):
        """Connect a button to a processor callback with input path logging."""
        if background:
            button.clicked.connect(
                lambda checked=False, name=action_name, cb=callback, control=button:
                self.run_background_processor_action(name, cb, control)
            )
        else:
            button.clicked.connect(
                lambda checked=False, name=action_name, cb=callback: self.run_processor_action(
                    name, cb
                )
            )

    def run_processor_action(self, action_name, callback):
        """Log current input selections before running a processor action."""
        self.log_selected_input_paths(action_name)
        self.save_input_state()
        callback()
        if action_name == "Fill DEM":
            self.update_l0_resolution_from_dem()

    def run_background_processor_action(self, action_name, callback, button):
        """Dispatch supported file jobs without moving QGIS objects off-thread."""
        self.log_selected_input_paths(action_name)
        self.save_input_state()
        controls = (button,)
        failed = lambda message: self._background_processor_failed(
            action_name, message
        )
        if action_name == "Fill DEM":
            started = self.morphology_tasks.start_fill(
                controls=controls, load=True, failed=failed
            )
        elif action_name == "Create Channel Network":
            started = self.morphology_tasks.start_hydrology(
                key="create-channel-network",
                controls=controls,
                load="channel_network",
                failed=failed,
            )
        elif action_name == "Flow Accumulation":
            started = self.morphology_tasks.start_hydrology(
                key="flow-accumulation",
                controls=controls,
                load="flow_accumulation",
                direction=False,
                include_channel=False,
                failed=failed,
            )
        elif action_name == "Flow Direction":
            started = self.morphology_tasks.start_hydrology(
                key="flow-direction",
                controls=controls,
                load="flow_direction",
                include_channel=False,
                failed=failed,
            )
        elif action_name in {"Slope", "Aspect"}:
            started = self.morphology_tasks.start_terrain(
                key=f"terrain-{action_name.lower()}",
                controls=controls,
                load=action_name.lower(),
                failed=failed,
            )
        elif action_name in {"Land Use", "Soil", "Hydrogeology"}:
            kind = {"Land Use": "lc", "Soil": "soil", "Hydrogeology": "geology"}[
                action_name
            ]
            started = self.morphology_tasks.start_categorical(
                kind,
                key=f"categorical-{kind}",
                controls=controls,
                load=True,
                failed=failed,
            )
        else:
            # Legacy callbacks remain on Qt's GUI thread until split into a
            # path-only worker and a main-thread completion callback.
            button.setEnabled(False)

            def run_on_main_thread():
                try:
                    callback()
                except Exception as error:
                    self._background_processor_failed(action_name, error)
                finally:
                    button.setEnabled(True)

            QtCore.QTimer.singleShot(0, run_on_main_thread)
            started = True
        if not started:
            QMessageBox.information(
                self, action_name, f"{action_name} is already running."
            )

    def _background_processor_failed(self, action_name, message):
        detail = str(message).split("\n", 1)[0]
        self.log_message(f"ERROR: {action_name}: {detail}")
        QMessageBox.critical(self, action_name, detail)

    def reset_geometry_processing(self):
        """Reset geometry outputs and refresh workflow UI state."""
        self.morphology_processor.resetGeometry()
        self.refresh_morphology_workflow_button_states()

    def handle_domain_definition_type(self, index):
        """Dispatch the selected domain-definition workflow."""
        if self._loading_input_state:
            return
        previous = self._domain_definition_mode
        previous_dem = self.checkBox_DEMdomain.isChecked()
        if int(index) == 2:
            self._domain_definition_mode = self.comboBox_domainDefinitionType.currentText()
            self.save_input_state()
            self.open_domain_delineator()
            return
        mode = "dem" if int(index) == 0 else "snapped"
        self.checkBox_DEMdomain.setChecked(mode == "dem")
        if self.open_domain_assignment(mode):
            self._domain_definition_mode = self.comboBox_domainDefinitionType.currentText()
            self.save_input_state()
            return
        self.checkBox_DEMdomain.setChecked(previous_dem)
        combo = self.comboBox_domainDefinitionType
        combo.blockSignals(True)
        combo.setCurrentIndex(combo.findText(previous) if previous else -1)
        combo.blockSignals(False)
        self.save_input_state()

    def open_domain_assignment(self, mode):
        """Open the non-interactive domain/gauge assignment workflow."""
        if standalone_is_active():
            QMessageBox.information(
                self,
                "Domain Assignment",
                "Domain preparation requires the QGIS plugin runtime.",
            )
            return False
        if not self.project_folder:
            QMessageBox.warning(self, "Domain Assignment", "Select a project folder first.")
            return False
        try:
            outlet_ids = self.selected_outlet_ids()
        except StationIdError as error:
            QMessageBox.warning(self, "Domain Assignment", str(error))
            return False
        if not outlet_ids:
            QMessageBox.warning(
                self,
                "Domain Assignment",
                "The pour-point layer does not contain any outlet features.",
            )
            return False
        try:
            from .Morphology.hydrology.discharge_dialog import (
                DischargeTableAssignmentDialog,
                DomainAndDischargeTableAssignmentDialog,
            )
            from .Morphology.watershed.domain_workflow import DomainWorkflow
            from .Morphology.watershed.domain_state import (
                DOMAIN_MODE_DEM_EXTENT,
                DOMAIN_MODE_SNAPPED,
            )
            from .nml_settings import sync_domain_settings

            layer = self.input_combo("pour_points").currentLayer()
            workflow = DomainWorkflow(
                self,
                self.morphology_processor,
                layer,
                self.selected_outlet_id_field(),
                outlet_ids,
            )
            definition = (
                DOMAIN_MODE_DEM_EXTENT if mode == "dem" else DOMAIN_MODE_SNAPPED
            )
            state = workflow.load_synced_state(definition, mode == "dem")
            dialog_class = (
                DischargeTableAssignmentDialog
                if mode == "dem"
                else DomainAndDischargeTableAssignmentDialog
            )
            dialog = dialog_class(
                outlet_ids,
                self,
                initial_records=state.get("outlets", {}),
            )
            execute = getattr(dialog, "exec", None) or dialog.exec_
            if execute() != QDialog.Accepted:
                return False
            assignments = dialog.selected_assignments()
            if mode == "dem":
                workflow.apply_dem_extent(assignments)
            else:
                workflow.apply_snapped_domains(assignments)
            sync_domain_settings(self.project_folder)
            self.morphology_processor.update_gauged_outlet_count()
            self.invalidate_meteo_morph_setup()
            return True
        except Exception as error:
            self.log_message(f"ERROR: Domain assignment failed: {error}")
            self.log_message(f"Traceback: {traceback.format_exc()}")
            QMessageBox.critical(self, "Domain Assignment", str(error))
            return False



    def open_domain_delineator(self, checked=False):
        """Open the per-outlet domain delineation dialog in QGIS."""
        if standalone_is_active():
            QMessageBox.information(
                self,
                "Domain Delineator",
                "Interactive domain delineation requires the QGIS plugin runtime.",
            )
            return
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Domain Delineator",
                "Select a project folder first.",
            )
            return
        try:
            outlet_ids = self.selected_outlet_ids()
        except StationIdError as error:
            QMessageBox.warning(self, "Domain Delineator", str(error))
            return
        if not outlet_ids:
            QMessageBox.warning(
                self,
                "Domain Delineator",
                "The pour-point layer does not contain any outlet features.",
            )
            return

        from .domain_delineator_dialog import DomainDelineatorDialog

        self.morphology_processor.load_project_state()
        layer = self.input_combo("pour_points").currentLayer()
        field = self.selected_outlet_id_field()
        controls = (
            self.comboBox_domainDefinitionType,
            self.pushButton_delineate,
        )
        started = self.morphology_tasks.start_domain_preflight(
            controls,
            lambda context: self._show_prepared_domain_delineator(
                context, layer, field, outlet_ids
            ),
            self._domain_delineator_preflight_failed,
        )
        if not started:
            QMessageBox.information(
                self,
                "Domain Delineator",
                "Domain Delineator preparation is already running.",
            )

    def _show_prepared_domain_delineator(
        self, context, pour_points_layer, outlet_id_field, outlet_ids
    ):
        from .domain_delineator_dialog import DomainDelineatorDialog

        try:
            dialog = DomainDelineatorDialog(
                self,
                self.morphology_processor,
                pour_points_layer,
                outlet_id_field,
                outlet_ids,
                prepared_context=context,
            )
            self._domain_delineator_dialog = dialog
            execute = getattr(dialog, "exec", None) or dialog.exec_
            execute()
        except Exception as error:
            self.log_message(f"ERROR: Domain delineator failed: {error}")
            self.log_message(f"Traceback: {traceback.format_exc()}")
            QMessageBox.critical(self, "Domain Delineator", str(error))
        finally:
            self._domain_delineator_dialog = None
            try:
                from .nml_settings import sync_domain_settings

                sync_domain_settings(self.project_folder)
            except Exception as error:
                self.log_message(
                    f"WARNING: Could not update domain namelist settings: {error}"
                )
            self.save_input_state()
            self.morphology_processor.update_gauged_outlet_count()

    def _domain_delineator_preflight_failed(self, message):
        message = str(message).split("\n", 1)[0]
        self.log_message(f"ERROR: Domain Delineator preparation failed: {message}")
        QMessageBox.critical(self, "Domain Delineator", message)

    def morphology_workflow_specs(self):
        """Return metadata for threaded morphology workflows."""
        return {
            "execute_all": {
                "button": "pushButton_executeAllMorphology",
                "fallback_button": "pushButton_executeAllMorphology",
                "label": "Execute All Processing",
                "action_name": "Execute All Processing",
                "method": "execute_all_processing",
                "thread_name": "mHM QGISExecuteAllThread",
                "completed_message": (
                    "Morphology preparation completed successfully."
                ),
                "failed_message": (
                    "Execute All Processing failed. Check the log for details."
                ),
            },
            "meteo_morph_setup": {
                "button": "pushButton_executeMeteoMorphSetup",
                "label": "Meteorology and Morphology Setup",
                "action_name": "Meteorology and Morphology Setup",
                "method": "execute_morph_setup_processing",
                "thread_name": "mHM QGISMeteoMorphSetupThread",
                "completed_message": (
                    "Meteorology and morphology setup completed successfully."
                ),
                "failed_message": (
                    "Meteorology and morphology setup failed. "
                    "Check the log for details."
                ),
            },
        }

    def _capture_workflow_button_default_styles(self):
        """Remember original button styles so saved states can be cleared."""
        for workflow_key in self.morphology_workflow_specs():
            button = self.morphology_workflow_button(workflow_key)
            if button is not None:
                self._workflow_button_default_styles[workflow_key] = (
                    button.styleSheet()
                )

    def morphology_workflow_button(self, workflow_key):
        """Return the button associated with a morphology workflow."""
        spec = self.morphology_workflow_specs().get(workflow_key, {})
        for button_name in (spec.get("button"), spec.get("fallback_button")):
            if not button_name:
                continue
            button = getattr(self, button_name, None)
            if button is not None:
                return button
        return None

    def running_morphology_workflow_key(self):
        """Return the first running morphology workflow key, if any."""
        if self.morphology_tasks.execute_all_active:
            return "execute_all"
        for workflow_key, thread in self._morphology_workflow_threads.items():
            if thread is not None and thread.isRunning():
                return workflow_key
        return None

    def start_execute_all_processing(self):
        """Run Execute All as a sequence of QGIS-managed file tasks."""
        self.morphology_tasks.start_execute_all()

    def start_meteo_morph_setup_processing(self):
        """Validate inputs and run meteo then morphology setup."""
        if self.running_morphology_workflow_key() is not None:
            self.start_morphology_workflow("meteo_morph_setup")
            return
        if not self.check_prerequisites():
            return

        try:
            precipitation, temperature, pet = self.selected_meteo_specs()
            inspections = inspect_meteo_inputs(
                precipitation,
                temperature,
                pet,
            )
            self._meteo_inspections = inspections
            grid = self.prepare_meteo_l2_grid(
                inspections["precipitation"]
            )
            self.set_meteo_l2_grid_metadata(grid["metadata"])
            target_grid = TargetGrid(
                lon=tuple(grid["lon"]),
                lat=tuple(grid["lat"]),
                header=dict(grid["header"]),
                crs=self.selected_meteo_crs() or None,
                sample_lon=grid["sample_lon"],
                sample_lat=grid["sample_lat"],
            )
            target_grid.validate()
            self._pending_meteo_run = MeteorologyRun(
                project_folder=Path(self.project_folder),
                precipitation=precipitation,
                temperature=temperature,
                pet=pet,
                target_grid=target_grid,
                grid_metadata=dict(grid["metadata"]),
            )
        except Exception as error:
            self.log_message(f"ERROR: Cannot start meteo setup: {error}")
            QMessageBox.warning(
                self,
                "Meteorology Setup",
                str(error),
            )
            return

        self.start_morphology_workflow("meteo_morph_setup")

    def start_morphology_workflow(self, workflow_key):
        """Start a named morphology workflow in a background worker thread."""
        if workflow_key == "execute_all":
            self.morphology_tasks.start_execute_all()
            return
        spec = self.morphology_workflow_specs().get(workflow_key)
        if spec is None:
            self.log_message(f"ERROR: Unknown morphology workflow: {workflow_key}")
            return

        if self.task_coordinator.resource_busy("morphology-processor"):
            QMessageBox.information(
                self,
                spec["label"],
                "Another morphology preprocessing task is currently running.",
            )
            return
        if not self.task_coordinator.has_capacity():
            QMessageBox.information(
                self,
                spec["label"],
                "All configured worker threads are currently busy.",
            )
            return

        running_key = self.running_morphology_workflow_key()
        if running_key is not None:
            running_spec = self.morphology_workflow_specs().get(
                running_key, {}
            )
            QMessageBox.information(
                self,
                spec["label"],
                f"{running_spec.get('label', 'A morphology workflow')} is already running.",
            )
            return

        if not self.check_prerequisites():
            return

        self.log_selected_input_paths(spec["action_name"])
        self.save_input_state()
        self.morphology_processor.load_processing_state()
        self.morphology_processor.mark_workflow_status(
            workflow_key,
            "running",
        )
        self.set_morphology_workflow_button_state(workflow_key, "running")
        if workflow_key == "meteo_morph_setup":
            self.set_meteo_setup_controls_enabled(False)

        thread = QtCore.QThread(self)
        thread.setObjectName(spec["thread_name"])
        worker = MeteorologyWorkflowWorker(
            workflow_key,
            spec["label"],
            meteorology_processor=self.meteorology_processor,
            meteorology_run=self._pending_meteo_run,
        )
        worker.moveToThread(thread)

        worker.log_message.connect(self.log_message)
        worker.finished.connect(self.finish_meteo_workflow_stage)
        worker.finished.connect(
            lambda key, ok, message, workflow_thread=thread: workflow_thread.quit()
        )
        worker.finished.connect(worker.deleteLater)
        thread.started.connect(worker.run)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(
            lambda key=workflow_key: self.clear_morphology_workflow_worker(key)
        )

        self._morphology_workflow_threads[workflow_key] = thread
        self._morphology_workflow_workers[workflow_key] = worker
        self.task_coordinator.start_external(
            f"workflow-{workflow_key}",
            spec["label"],
            resource="morphology-processor",
        )
        worker.log_message.connect(
            lambda message, key=workflow_key: self.task_coordinator.append_log(
                f"workflow-{key}", message
            )
        )
        thread.start()

    def finish_meteo_workflow_stage(self, workflow_key, ok, message):
        """Run QGIS-dependent morphology on the GUI thread after meteorology."""
        if not ok:
            self.finish_morphology_workflow(workflow_key, False, message)
            return

        self.morphology_processor.load_processing_state()
        try:
            ok = bool(self.morphology_processor.execute_morph_setup_processing(
                show_error_dialog=False,
                workflow_key=workflow_key,
            ))
            self.finish_morphology_workflow(workflow_key, ok, "")
        except Exception as error:
            self.log_message(
                f"ERROR: Morphology setup failed with exception: {error}\n"
                f"{traceback.format_exc()}"
            )
            self.finish_morphology_workflow(
                workflow_key,
                False,
                str(error),
            )

    def finish_morphology_workflow(self, workflow_key, ok, message):
        """Update UI and persisted workflow status after a workflow finishes."""
        spec = self.morphology_workflow_specs().get(workflow_key)
        if spec is None:
            return
        if workflow_key == "meteo_morph_setup":
            self.set_meteo_setup_controls_enabled(True)

        if ok:
            self.morphology_processor.mark_workflow_status(
                workflow_key,
                "completed",
                spec["completed_message"],
            )
            self.set_morphology_workflow_button_state(workflow_key, "completed")
            self.log_message(spec["completed_message"])
            self.task_coordinator.finish_external(
                f"workflow-{workflow_key}", True, spec["completed_message"]
            )
            return

        workflow_message = self.morphology_processor.workflow_status(
            workflow_key
        ).get("message")
        message = (
            message
            or workflow_message
            or spec["failed_message"]
        )
        self.morphology_processor.mark_workflow_status(
            workflow_key,
            "failed",
            message,
        )
        self.set_morphology_workflow_button_state(workflow_key, "failed")
        self.log_message(message)
        self.task_coordinator.finish_external(
            f"workflow-{workflow_key}", False, message
        )
        QMessageBox.critical(
            self,
            spec["label"],
            message,
        )

    def clear_morphology_workflow_worker(self, workflow_key):
        """Drop worker references after a morphology workflow thread stops."""
        self._morphology_workflow_threads.pop(workflow_key, None)
        self._morphology_workflow_workers.pop(workflow_key, None)
        if workflow_key == "meteo_morph_setup":
            self._pending_meteo_run = None

    def invalidate_meteo_morph_setup(self):
        """Clear saved completion when any meteo setup input changes."""
        workflow_key = "meteo_morph_setup"
        if (
                self._loading_input_state
                or not self.project_folder
                or self.running_morphology_workflow_key() == workflow_key):
            return
        self.morphology_processor.load_processing_state()
        workflows = self.morphology_processor.processing_state.setdefault(
            "workflows", {}
        )
        if workflow_key in workflows:
            workflows.pop(workflow_key)
            self.morphology_processor.save_processing_state()
        self.set_morphology_workflow_button_state(workflow_key, "")

    def set_meteo_setup_controls_enabled(self, enabled):
        """Freeze all run-affecting dialog controls during the combined run."""
        self.pushButton_BrowseProjectFolder.setEnabled(enabled)
        self.tabWidget.setEnabled(enabled)
        self.stackedWidget.setEnabled(enabled)

    def closeEvent(self, event):
        """Prevent closing the dialog while a morphology workflow is running."""
        running_key = self.running_morphology_workflow_key()
        if running_key is not None or self.task_coordinator.is_busy():
            spec = self.morphology_workflow_specs().get(running_key, {})
            QMessageBox.warning(
                self,
                spec.get("label", "Morphology Processing"),
                (
                    f"{spec.get('label', 'Morphology processing')} is still running. "
                    "Wait for it to finish before closing mHM QGIS."
                ),
            )
            event.ignore()
            return
        super().closeEvent(event)

    def set_execute_all_button_state(self, status):
        """Reflect execute-all state on the toolbar-style button."""
        self.set_morphology_workflow_button_state("execute_all", status)

    def set_morphology_workflow_button_state(self, workflow_key, status):
        """Reflect workflow state on its toolbar-style button."""
        button = self.morphology_workflow_button(workflow_key)
        if button is None:
            return

        if status == "running":
            button.setEnabled(False)
            button.setStyleSheet(
                "QPushButton {"
                "text-align: left;"
                "background-color: #f6c453;"
                "border: 1px solid #a66f00;"
                "border-radius: 3px;"
                "}"
            )
            return

        button.setEnabled(True)
        if status == "completed":
            button.setStyleSheet(
                "QPushButton {"
                "text-align: left;"
                "background-color: #2e7d32;"
                "border: 1px solid #1b5e20;"
                "border-radius: 3px;"
                "}"
            )
        elif status == "failed":
            button.setStyleSheet(
                "QPushButton {"
                "text-align: left;"
                "background-color: #c62828;"
                "border: 1px solid #8e0000;"
                "border-radius: 3px;"
                "}"
            )
        else:
            button.setStyleSheet(
                self._workflow_button_default_styles.get(workflow_key, "")
            )

    def refresh_execute_all_button_state(self):
        """Restore execute-all button styling from project processing state."""
        self.refresh_morphology_workflow_button_state("execute_all")

    def refresh_morphology_workflow_button_states(self):
        """Restore all morphology workflow button styles from project state."""
        for workflow_key in self.morphology_workflow_specs():
            self.refresh_morphology_workflow_button_state(workflow_key)

    def refresh_morphology_workflow_button_state(self, workflow_key):
        """Restore one morphology workflow button style from project state."""
        if self.running_morphology_workflow_key() == workflow_key:
            return

        workflow = self.morphology_processor.workflow_status(workflow_key)
        if workflow.get("status") == "completed":
            if (
                    workflow_key == "meteo_morph_setup"
                    and self.project_folder
                    and not os.path.isfile(
                        os.path.join(
                            data_folder(self.project_folder),
                            "latlon.nc",
                        )
                    )):
                self.invalidate_meteo_morph_setup()
                return
            self.set_morphology_workflow_button_state(workflow_key, "completed")
        elif workflow.get("status") == "failed":
            self.set_morphology_workflow_button_state(workflow_key, "failed")
        else:
            self.set_morphology_workflow_button_state(workflow_key, "")

    def project_terminal_dialog(self):
        """Return the persistent terminal dialog for this plugin dialog."""
        if self._terminal_dialog is None:
            self._terminal_dialog = ProjectTerminalDialog(self)
        return self._terminal_dialog

    def open_thread_display(self, checked=False):
        """Show the persistent task monitor for plugin preprocessing."""
        if self._thread_display_dialog is None:
            self._thread_display_dialog = ThreadDisplayDialog(
                self.task_coordinator, self
            )
        self._thread_display_dialog.show()
        self._thread_display_dialog.raise_()
        self._thread_display_dialog.activateWindow()
        return self._thread_display_dialog

    def open_project_terminal(self):
        """Open the persistent terminal in the plugin-owned workspace."""
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Project Folder Required",
                "Select a project folder before opening the terminal.",
            )
            return None
        terminal = self.project_terminal_dialog()
        terminal.show_for_directory(workspace_folder(self.project_folder))
        return terminal

    def categorical_type_combo(self, kind):
        """Return the land-cover, soil, or geology input-type combo box."""
        name = {
            "lc": "comboBox_landUseInputType",
            "soil": "comboBox_soil_inputType",
            "geology": "comboBox_geology_inputType",
            "lai": "comboBox_lai_inputType",
        }[kind]
        combo = getattr(self, name, None)
        if combo is None and kind == "lc":
            combo = self.comboBox_landUseInputType
        return combo

    def categorical_input_mode(self, kind):
        """Return the selected categorical input mode."""
        return self.categorical_type_combo(kind).currentText().strip()

    def categorical_lookup_config(self, kind):
        """Return the accepted lookup-table selection for one data type."""
        config = self._categorical_lookup_configs.get(kind)
        return dict(config) if config else None

    def categorical_source_config(self, kind):
        """Return the dialog-owned source for a categorical workflow."""
        if self.categorical_input_mode(kind).strip().lower() == "mhm ready":
            config = self._categorical_ready_configs.get(kind)
        else:
            config = self._categorical_lookup_configs.get(kind)
        return dict(config) if isinstance(config, dict) else None

    def lai_netcdf_config(self):
        return dict(self._lai_input_config)

    def handle_categorical_type(self, kind, text):
        """Open the lookup dialog immediately when lookup mode is selected."""
        text = str(text or "").strip()
        previous = self._categorical_modes.get(kind, "")
        if self._loading_input_state:
            self._categorical_modes[kind] = text
            return
        normalized = text.lower()
        if kind == "lai":
            self.handle_lai_input_type(text, previous)
            return
        if kind == "lc":
            self.handle_land_use_input_type(text, previous)
            return
        if kind == "soil" and "multi-horizon" in normalized:
            self.handle_multi_horizon_soil_input(text, previous)
            return
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Project Folder Required",
                "Select a project folder before configuring a lookup table.",
            )
            self._restore_categorical_mode(kind, previous)
            return

        if normalized == "mhm ready":
            dialog = MhmReadyInputDialog(
                self.project_folder,
                kind,
                self,
                initial=self._categorical_ready_configs.get(kind),
            )
        elif "lookup table" in normalized:
            dialog = SingleLayerInputDialog(
                self.project_folder,
                kind,
                self,
                initial=self._categorical_lookup_configs.get(kind),
            )
        else:
            self._categorical_modes[kind] = text
            self.save_input_state()
            self.invalidate_meteo_morph_setup()
            return
        execute = getattr(dialog, "exec", None) or dialog.exec_
        if execute() != QDialog.Accepted or dialog.selected_config() is None:
            self._restore_categorical_mode(kind, previous)
            return

        config = dialog.selected_config()
        if normalized == "mhm ready":
            self._categorical_ready_configs[kind] = {
                "input_path": config.input_path,
                "classdefinition_path": config.classdefinition_path,
            }
            self._categorical_lookup_configs.pop(kind, None)
        else:
            self._categorical_lookup_configs[kind] = {
                "input_path": config.input_path,
                "lookup_table": config.lookup_table,
                "mapping_field": config.mapping_field,
                "class_field": config.class_field,
            }
            self._categorical_ready_configs.pop(kind, None)
        self._categorical_modes[kind] = text
        if kind == "soil":
            self._advanced_inputs.pop("soil", None)
            self._save_standard_soil_nml_input(
                "mhm_ready" if normalized == "mhm ready" else "single_categorical"
            )
        self.save_input_state()
        self.invalidate_meteo_morph_setup()

    def handle_lai_input_type(self, text, previous=""):
        """Collect LAI NetCDF or future categorical input configuration."""
        if not self.project_folder:
            QMessageBox.warning(
                self, "Project Folder Required", "Select a project folder first."
            )
            self._restore_categorical_mode("lai", previous)
            return
        normalized = str(text or "").strip().lower()
        if "netcdf" in normalized:
            dialog = LaiNetcdfInputDialog(
                self.project_folder,
                self,
                initial=self._lai_input_config,
            )
        else:
            dialog = SingleLayerInputDialog(
                self.project_folder,
                "lai",
                self,
                initial=self._categorical_lookup_configs.get("lai"),
            )
        execute = getattr(dialog, "exec", None) or dialog.exec_
        if execute() != QDialog.Accepted or dialog.selected_config() is None:
            self._restore_categorical_mode("lai", previous)
            return
        config = dialog.selected_config()
        if "netcdf" in normalized:
            self._lai_input_config = {
                "input_path": config.input_path,
                "input_resolution": config.input_resolution,
                "target_timestep": config.target_timestep,
            }
            self._categorical_lookup_configs.pop("lai", None)
            from .lai_temporal import lai_time_step
            from .nml_settings import update_section

            update_section(
                self.project_folder,
                "lai",
                {
                    "mode": "netcdf",
                    "source_path": config.input_path,
                    "source_variable": "",
                    "input_resolution": config.input_resolution,
                    "target_timestep": config.target_timestep,
                    "time_step": lai_time_step(config.target_timestep),
                    "output_path": "data/master/lai/lai.nc",
                    "variable": "lai",
                },
            )
        else:
            self._categorical_lookup_configs["lai"] = {
                "input_path": config.input_path,
                "lookup_table": config.lookup_table,
                "mapping_field": config.mapping_field,
                "class_field": config.class_field,
            }
            self._lai_input_config = {}
            from .nml_settings import update_section

            update_section(
                self.project_folder,
                "lai",
                {
                    "mode": "single_categorical",
                    **self._categorical_lookup_configs["lai"],
                    "processing_status": "not_implemented",
                },
            )
        self._categorical_modes["lai"] = str(text)
        self.save_input_state()
        self.invalidate_meteo_morph_setup()

    def handle_land_use_input_type(self, text, previous=""):
        """Collect the selected ready, single, or historical land-use input."""
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Project Folder Required",
                "Select a project folder before configuring land use.",
            )
            self._restore_categorical_mode("lc", previous)
            return
        normalized = str(text).strip().lower()
        if normalized == "mhm ready":
            dialog = MhmReadyInputDialog(
                self.project_folder,
                "lc",
                self,
                initial=self._categorical_ready_configs.get("lc"),
            )
            execute = getattr(dialog, "exec", None) or dialog.exec_
            if execute() != QDialog.Accepted or dialog.selected_config() is None:
                self._restore_categorical_mode("lc", previous)
                return
            config = dialog.selected_config()
            self._categorical_ready_configs["lc"] = {
                "input_path": config.input_path,
                "classdefinition_path": "",
            }
            self._categorical_lookup_configs.pop("lc", None)
            self._advanced_inputs.pop("land_cover", None)
            self._land_cover_ready_source = config.input_path
        elif "single categorical" in normalized:
            dialog = SingleLayerInputDialog(
                self.project_folder,
                "lc",
                self,
                initial=self._categorical_lookup_configs.get("lc"),
            )
            execute = getattr(dialog, "exec", None) or dialog.exec_
            if execute() != QDialog.Accepted or dialog.selected_config() is None:
                self._restore_categorical_mode("lc", previous)
                return
            config = dialog.selected_config()
            self._categorical_lookup_configs["lc"] = {
                "input_path": config.input_path,
                "lookup_table": config.lookup_table,
                "mapping_field": config.mapping_field,
                "class_field": config.class_field,
            }
            self._categorical_ready_configs.pop("lc", None)
            self._advanced_inputs.pop("land_cover", None)
            self._land_cover_ready_source = ""
        else:
            from .advanced_input_dialogs import HistoricalLandUseDialog

            initial = self._advanced_inputs.get("land_cover")
            dialog = HistoricalLandUseDialog(
                self.project_folder,
                self,
                initial=initial,
                all_time="single layer" in str(text).lower(),
            )
            execute = getattr(dialog, "exec", None) or dialog.exec_
            if execute() != QDialog.Accepted:
                self._restore_categorical_mode("lc", previous)
                return
            value = dialog.selected_input()
            self._advanced_inputs["land_cover"] = value
            self._land_cover_ready_source = ""
            self._save_land_cover_nml_input(value)
        self._categorical_modes["lc"] = str(text)
        self.save_input_state()
        self.invalidate_meteo_morph_setup()
        self.update_morphology_date_control()

    def handle_multi_horizon_soil_input(self, text, previous=""):
        """Collect and persist the multi-horizon soil source configuration."""
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Project Folder Required",
                "Select a project folder before configuring soil inputs.",
            )
            self._restore_categorical_mode("soil", previous)
            return
        from .advanced_input_dialogs import MultiHorizonSoilDialog

        dialog = MultiHorizonSoilDialog(
            self.project_folder,
            self,
            initial=self._advanced_inputs.get("soil"),
        )
        execute = getattr(dialog, "exec", None) or dialog.exec_
        if execute() != QDialog.Accepted:
            self._restore_categorical_mode("soil", previous)
            return
        value = dialog.selected_input()
        self._advanced_inputs["soil"] = value
        self._categorical_modes["soil"] = str(text)
        self._save_soil_nml_input(value)
        self.save_input_state()
        self.invalidate_meteo_morph_setup()

    def uses_advanced_categorical_input(self, kind):
        """Return whether processing should use the advanced input pipeline."""
        text = self.categorical_input_mode(kind).lower()
        return (
            (kind == "lc" and "historical" in text)
            or (kind == "soil" and "multi-horizon" in text)
        )

    def process_advanced_categorical_input(self, kind):
        """Prepare configured historical land use or multi-horizon soil."""
        if not self.check_prerequisites():
            return False
        if not self.morphology_processor._ensure_filled_dem(
            self.morphology_processor.fill_dem
        ):
            return False
        version = self.comboBox_mHMversion.currentText().strip()
        try:
            if kind == "lc" and self.categorical_input_mode("lc").lower() == "mhm ready":
                from .advanced_input_processing import configure_ready_land_cover

                if not self._land_cover_ready_source:
                    raise ValueError("Select an mHM-ready land-cover file first.")
                outputs = (
                    configure_ready_land_cover(
                        self.project_folder,
                        self._land_cover_ready_source,
                        version,
                    ),
                )
            elif kind == "lc":
                from .advanced_input_processing import process_land_cover_input

                value = self._land_use_input_value(
                    self._advanced_inputs.get("land_cover")
                )
                outputs = process_land_cover_input(
                    self.project_folder,
                    version,
                    value,
                    self.morphology_processor.filled_dem_path,
                    log=self.log_message,
                )
            else:
                from .advanced_input_processing import process_soil_input

                value = self._soil_input_value(self._advanced_inputs.get("soil"))
                outputs = process_soil_input(
                    self.project_folder,
                    version,
                    value,
                    self.morphology_processor.filled_dem_path,
                    log=self.log_message,
                )
            for output in outputs:
                self.morphology_processor.mark_output_prepared(
                    str(output), name=output.name, loaded=False
                )
            self.log_message(
                f"Advanced {'land-cover' if kind == 'lc' else 'soil'} data prepared."
            )
            self.update_morphology_date_control()
            return True
        except Exception as error:
            self.log_message(f"ERROR preparing advanced {kind} data: {error}")
            QMessageBox.critical(self, "Morphology Input Error", str(error))
            return False

    def _save_land_cover_nml_input(self, value):
        from .nml_settings import update_section

        update_section(
            self.project_folder,
            "land_cover",
            {
                "mode": "historical" if len(value.periods) > 1 else "single",
                "variable": "land_cover",
                "lookup_table": str(value.lookup_table),
                "mapping_field": value.mapping_field,
                "class_field": value.class_field,
                "scenes": [
                    {
                        "start_year": item.start_year,
                        "end_year": item.end_year,
                        "source_path": str(item.file_path),
                    }
                    for item in value.periods
                ],
            },
        )

    def _save_soil_nml_input(self, value):
        from .nml_settings import update_section

        data = value.as_dict()
        data.update(
            {
                "mode": "multi_horizon",
                "soil_db_mode": (
                    0 if self.comboBox_mHMversion.currentText().startswith("5") else 1
                ),
                "variable": "soil_class",
                "source_bulk_density_unit": value.bulk_density_unit,
                "bulk_density_unit": "g/cm3",
                "composition_normalization": "component_sum_percent",
            }
        )
        update_section(self.project_folder, "soil", data)

    def _save_standard_soil_nml_input(
        self, mode, output_path=None, classdefinition_path=None
    ):
        if not self.project_folder:
            return
        from .nml_settings import relative_workspace_path, update_section

        output = (
            relative_workspace_path(self.project_folder, output_path)
            if output_path
            else "data/master/static/morph/soil_class.asc"
        )
        definition = (
            relative_workspace_path(self.project_folder, classdefinition_path)
            if classdefinition_path
            else "data/master/static/morph/soil_classdefinition.txt"
        )
        update_section(
            self.project_folder,
            "soil",
            {
                "mode": mode,
                "soil_db_mode": 0,
                "output_path": output,
                "variable": "soil_class",
                "classdefinition_path": definition,
                "horizons": [],
            },
        )

    def record_standard_soil_output(self, output_path=None, classdefinition_path=None):
        """Update the namelist handoff after standard soil preparation."""
        self._save_standard_soil_nml_input(
            "mhm_ready" if self.categorical_input_mode("soil").lower() == "mhm ready"
            else "single_categorical",
            output_path,
            classdefinition_path,
        )

    def refresh_advanced_nml_settings(self):
        """Refresh version-dependent settings after the mHM version changes."""
        value = self._advanced_inputs.get("soil")
        if value:
            self._save_soil_nml_input(self._soil_input_value(value))
        land_cover = self._advanced_inputs.get("land_cover")
        if land_cover:
            self._save_land_cover_nml_input(
                self._land_use_input_value(land_cover)
            )
        self.invalidate_meteo_morph_setup()

    @staticmethod
    def _land_use_input_value(value):
        from .advanced_input_manifests import LandUseInput, LandUsePeriod

        if isinstance(value, LandUseInput):
            return value
        if not isinstance(value, dict):
            raise ValueError("Configure the land-use layers first.")
        return LandUseInput(
            tuple(
                LandUsePeriod(
                    int(item["start_year"]),
                    int(item["end_year"]),
                    Path(item["file_path"]),
                )
                for item in value.get("periods", [])
            ),
            Path(value.get("lookup_table", "")),
            str(value.get("mapping_field", "")),
            str(value.get("class_field", "")),
        )

    @staticmethod
    def _soil_input_value(value):
        from .advanced_input_manifests import SoilHorizon, SoilInput

        if isinstance(value, SoilInput):
            return value
        if not isinstance(value, dict):
            raise ValueError("Configure the multi-horizon soil layers first.")
        return SoilInput(
            tuple(
                SoilHorizon(
                    int(item["horizon"]),
                    float(item["upper_depth"]),
                    float(item["lower_depth"]),
                    Path(item["clay_layer"]),
                    Path(item["sand_layer"]),
                    Path(item["silt_layer"]),
                    Path(item["bulk_density_layer"]),
                )
                for item in value.get("horizons", [])
            ),
            str(value.get("bulk_density_unit", "")),
        )

    def _restore_categorical_mode(self, kind, text):
        combo = self.categorical_type_combo(kind)
        combo.blockSignals(True)
        index = combo.findText(text) if text else -1
        combo.setCurrentIndex(index)
        combo.blockSignals(False)

    def restore_lai_input_type(self, lai_input_type):
        """Restore the selected LAI input type by text when possible."""
        combo_box = self.comboBox_lai_inputType
        if combo_box is None or not lai_input_type:
            return
        index = combo_box.findText(lai_input_type)
        if index >= 0:
            combo_box.setCurrentIndex(index)

    def input_layer_widgets(self):
        """Return persistent layer input widgets and their state keys."""
        widgets = []
        for key, kind in (
            ("dem", "dem"),
            ("pour_points", "pour_points"),
            ("soil", "soil"),
            ("land_cover", "land_cover"),
            ("geology", "geology"),
            ("lai_class", "lai"),
        ):
            widget = self.input_combo(kind)
            if widget is not None:
                widgets.append((key, widget))

        return widgets

    def input_text_widgets(self):
        """Return persistent text/path input widgets and their state keys."""
        widgets = []

        for key, widget_name in (
            ("lai_file", "lineEdit_lai_file"),
        ):
            widget = getattr(self, widget_name, None)
            if widget is not None:
                widgets.append((key, widget))
        return widgets

    def input_state_path(self):
        """Return the project-local input state file path."""
        if not self.project_folder:
            return None
        return os.path.join(
            workspace_folder(self.project_folder),
            self.input_state_filename,
        )

    def save_input_state(self):
        """Save selected inputs to a JSON file in the project folder."""
        if (
            not self.project_folder
            or self._loading_input_state
            or self._suspend_input_state_saves
        ):
            return

        state_path = self.input_state_path()
        if not state_path:
            return

        existing_state = self.read_existing_input_state(state_path)
        existing_layers = existing_state.get("layers", {})
        if not isinstance(existing_layers, dict):
            existing_layers = {}

        layers = {}
        for key, combo_box in self.input_layer_widgets():
            item_data = combo_box.currentData()
            if isinstance(item_data, dict) and item_data.get("origin") == "file":
                source = item_data.get("path", "")
                suffix = os.path.splitext(source)[1].lower()
                layers[key] = {
                    "name": os.path.basename(source),
                    "source": source,
                    "type": "vector" if suffix == ".shp" else "raster",
                    "origin": "file",
                    "manual": bool(item_data.get("manual")),
                }
                continue
            layer = combo_box.currentLayer()
            if not layer:
                if self._preserve_missing_layer_state and key in existing_layers:
                    layers[key] = existing_layers.get(key)
                else:
                    layers[key] = None
                continue

            layer_type = (
                "raster"
                if isinstance(layer, QgsRasterLayer)
                or layer.type() == QgsMapLayer.RasterLayer
                else "vector"
            )
            layers[key] = {
                "name": layer.name(),
                "source": layer.source(),
                "type": layer_type,
                "origin": "qgis",
            }

        text_inputs = {key: widget.text() for key, widget in self.input_text_widgets()}
        grid_resolutions = {
            "l1_resolution": self.current_l1_resolution(),
            "l11_resolution": self.current_l11_resolution(),
        }

        crs = self.get_crs()
        state = {
            "version": 4,
            "mhm_version": self.comboBox_mHMversion.currentText().strip(),
            "layers": layers,
            "text_inputs": text_inputs,
            "grid_resolutions": grid_resolutions,
            "grid_configuration": self.grid_configuration_snapshot(),
            "lai_input_type": self.comboBox_lai_inputType.currentText(),
            "folder_search": self.checkBox_enableFolderSearch.isChecked(),
            "categorical_types": {
                kind: self.categorical_input_mode(kind)
                for kind in ("lc", "soil", "geology", "lai")
            },
            "categorical_lookups": self._serialized_categorical_lookups(),
            "categorical_ready": self._serialized_path_configs(
                self._categorical_ready_configs
            ),
            "lai_netcdf": self._serialized_path_config(self._lai_input_config),
            "thread_count": self.task_coordinator.max_threads,
            "advanced_inputs": {
                kind: value.as_dict() if hasattr(value, "as_dict") else value
                for kind, value in self._advanced_inputs.items()
            },
            "land_cover_ready_source": self._land_cover_ready_source,
            "meteo_inputs": self.serialized_meteo_inputs(),
            "pour_point_outlet_id_field": self.selected_outlet_id_field(),
            "dem_domain": self.checkBox_DEMdomain.isChecked(),
            "domain_definition_type": self.comboBox_domainDefinitionType.currentText().strip(),
            "crs_authid": crs.authid() if crs and crs.isValid() else "",
            "project_layout": {
                "data_folder": data_folder(self.project_folder),
                "data_raw_folder": data_raw_folder(self.project_folder),
                "z_temp_folder": z_temp_folder(self.project_folder),
                "geometry_folder": geometry_folder(self.project_folder),
                "output_folder": output_folder(self.project_folder),
                "restart_folder": restart_folder(self.project_folder),
            },
        }

        try:
            os.makedirs(os.path.dirname(state_path), exist_ok=True)
            with open(state_path, "w", encoding="utf-8") as state_file:
                json.dump(state, state_file, indent=2, sort_keys=True)
        except Exception as e:
            self.log_message(f"WARNING: Could not save input state: {e}")

    def _serialized_categorical_lookups(self):
        configs = {}
        for kind, config in self._categorical_lookup_configs.items():
            saved = dict(config)
            for key in ("input_path", "lookup_table"):
                saved[key] = self._portable_input_path(saved.get(key))
            configs[kind] = saved
        return configs

    def _serialized_path_configs(self, configs):
        return {
            kind: self._serialized_path_config(config)
            for kind, config in configs.items()
            if isinstance(config, dict)
        }

    def _serialized_path_config(self, config):
        if not isinstance(config, dict):
            return {}
        saved = dict(config)
        for key in ("input_path", "classdefinition_path"):
            if key in saved:
                saved[key] = self._portable_input_path(saved.get(key))
        return saved

    def _portable_input_path(self, path):
        if not path or not self.project_folder:
            return str(path or "")
        try:
            relative = os.path.relpath(path, self.project_folder)
            if relative != ".." and not relative.startswith(f"..{os.sep}"):
                return relative.replace("\\", "/")
        except ValueError:
            pass
        return str(path)

    def serialized_meteo_inputs(self):
        """Return project-portable meteorology folder and source selections."""
        state = {
            "l2_multiplier": int(self.spinBox_L2ResolutionMultiplier.value()),
        }
        for kind, combo, source_combo in self.meteo_input_widgets():
            path = self.selected_folder_path(combo)
            if path and self.project_folder:
                try:
                    relative = os.path.relpath(path, self.project_folder)
                    if relative != ".." and not relative.startswith(f"..{os.sep}"):
                        path = relative.replace("\\", "/")
                except ValueError:
                    pass
            state[kind] = {
                "folder": path,
                "source": source_combo.currentText().strip(),
            }
        return state

    def read_existing_input_state(self, state_path):
        """Return existing input state for preservation during teardown."""
        if not state_path or not os.path.exists(state_path):
            return {}
        try:
            with open(state_path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if isinstance(state, dict):
                return state
        except Exception:
            pass
        return {}

    def load_input_state(self):
        """Load saved input selections from the project folder."""
        state_path = self.input_state_path()
        self._loading_input_state = True
        try:
            self._clear_input_state_selections()
            if not state_path or not os.path.exists(state_path):
                self.refresh_input_sources()
                self.refresh_meteo_folder_sources()
                self.log_message("No saved input state found for this project.")
                return
            try:
                with open(state_path, "r", encoding="utf-8") as state_file:
                    state = json.load(state_file)
            except Exception as error:
                self.refresh_input_sources()
                self.log_message(f"WARNING: Could not read input state: {error}")
                return

            self.restore_mhm_version(state.get("mhm_version", ""))
            self.restore_grid_resolution_preferences(
                state.get("grid_resolutions", {}))
            self.restore_text_inputs(state.get("text_inputs", {}))
            self.restore_lai_input_type(state.get("lai_input_type", ""))
            self.restore_input_crs(state.get("crs_authid", ""))
            self.restore_categorical_inputs(state)
            self.restore_advanced_inputs(state)
            self.refresh_input_sources()
            self.refresh_meteo_folder_sources()
            self.restore_input_layers(state.get("layers", {}))
            self.restore_domain_inputs(state)
            self.restore_meteo_inputs(state.get("meteo_inputs", {}))
            self.refresh_grid_resolution_controls()
        finally:
            self._loading_input_state = False

        for kind, _, _ in self.meteo_input_widgets():
            self.inspect_meteo_selection(kind, show_errors=False)
        self.log_message(f"Input state loaded: {state_path}")

    def _clear_input_state_selections(self):
        """Prevent selections from leaking between project folders."""
        for _, combo in self.input_layer_widgets():
            combo.setCurrentIndex(-1)
        for kind in ("lc", "soil", "geology", "lai"):
            combo = self.categorical_type_combo(kind)
            combo.blockSignals(True)
            combo.setCurrentIndex(-1)
            combo.blockSignals(False)
        for _, widget in self.input_text_widgets():
            widget.clear()
        self.comboBox_pourPointOutletID.blockSignals(True)
        self.comboBox_pourPointOutletID.clear()
        self.comboBox_pourPointOutletID.addItem("")
        self.comboBox_pourPointOutletID.blockSignals(False)
        self.checkBox_DEMdomain.blockSignals(True)
        self.checkBox_DEMdomain.setChecked(False)
        self.checkBox_DEMdomain.blockSignals(False)
        self.comboBox_domainDefinitionType.blockSignals(True)
        self.comboBox_domainDefinitionType.setCurrentIndex(-1)
        self.comboBox_domainDefinitionType.blockSignals(False)
        for _, folder_combo, source_combo in self.meteo_input_widgets():
            folder_combo.blockSignals(True)
            folder_combo.setCurrentIndex(-1)
            folder_combo.blockSignals(False)
            source_combo.blockSignals(True)
            source_combo.setCurrentIndex(-1)
            source_combo.blockSignals(False)
        self.spinBox_L2ResolutionMultiplier.blockSignals(True)
        self.spinBox_L2ResolutionMultiplier.setValue(1)
        self.spinBox_L2ResolutionMultiplier.blockSignals(False)
        self._meteo_inspections = {}
        self.clear_precipitation_resolution_labels()
        self.checkBox_enableFolderSearch.blockSignals(True)
        self.checkBox_enableFolderSearch.setChecked(False)
        self.checkBox_enableFolderSearch.blockSignals(False)
        self._categorical_lookup_configs = {}
        self._categorical_ready_configs = {}
        self._lai_input_config = {}
        self._categorical_modes = {}
        self._advanced_inputs = {}
        self._land_cover_ready_source = ""
        self._domain_definition_mode = ""
        self._preferred_l1_resolution = None
        self._preferred_l11_resolution = None
        self.task_coordinator.set_max_threads(2)

    def restore_domain_inputs(self, state):
        """Restore the selected outlet ID field and DEM-domain flag."""
        self.populate_pour_point_outlet_fields(
            preferred=str(state.get("pour_point_outlet_id_field", "") or "")
        )
        self.checkBox_DEMdomain.blockSignals(True)
        self.checkBox_DEMdomain.setChecked(bool(state.get("dem_domain", False)))
        self.checkBox_DEMdomain.blockSignals(False)
        combo = self.comboBox_domainDefinitionType
        if combo is not None:
            text = str(state.get("domain_definition_type", "") or "")
            combo.blockSignals(True)
            combo.setCurrentIndex(combo.findText(text) if text else -1)
            combo.blockSignals(False)
            self._domain_definition_mode = text

    def restore_meteo_inputs(self, state):
        """Restore the three meteo folders, source modes, and L2 multiplier."""
        if not isinstance(state, dict):
            return
        try:
            multiplier = max(1, int(state.get("l2_multiplier", 1)))
        except (TypeError, ValueError):
            multiplier = 1
        self.spinBox_L2ResolutionMultiplier.setValue(multiplier)

        for kind, folder_combo, source_combo in self.meteo_input_widgets():
            saved = state.get(kind, {})
            if not isinstance(saved, dict):
                continue
            folder = str(saved.get("folder", "") or "")
            if folder and not os.path.isabs(folder):
                folder = os.path.abspath(os.path.join(self.project_folder, folder))
            index = self._folder_combo_index(folder_combo, folder)
            if index < 0 and folder and os.path.isdir(folder):
                self._add_folder_combo_item(folder_combo, folder)
                index = folder_combo.count() - 1
            folder_combo.setCurrentIndex(index if index >= 0 else 0)

            source = str(saved.get("source", "") or "")
            source_index = source_combo.findText(source)
            if source_index < 0:
                normalized = source.lower().replace("_", " ")
                for candidate in range(source_combo.count()):
                    if source_combo.itemText(candidate).lower().replace("_", " ") == normalized:
                        source_index = candidate
                        break
            source_combo.setCurrentIndex(source_index)

    def restore_categorical_inputs(self, state):
        """Restore folder search, categorical modes, and lookup selections."""
        self.checkBox_enableFolderSearch.blockSignals(True)
        self.checkBox_enableFolderSearch.setChecked(bool(state.get("folder_search")))
        self.checkBox_enableFolderSearch.blockSignals(False)

        configs = {}
        for kind, config in state.get("categorical_lookups", {}).items():
            if kind not in {"lc", "soil", "geology", "lai"} or not isinstance(config, dict):
                continue
            restored = dict(config)
            for key in ("input_path", "lookup_table"):
                path = restored.get(key)
                if path and not os.path.isabs(path):
                    restored[key] = os.path.abspath(
                        os.path.join(self.project_folder, path)
                    )
            configs[kind] = restored
        self._categorical_lookup_configs = configs

        ready = {}
        for kind, config in state.get("categorical_ready", {}).items():
            if kind not in {"lc", "soil", "geology"} or not isinstance(config, dict):
                continue
            restored = dict(config)
            for key in ("input_path", "classdefinition_path"):
                path = restored.get(key)
                if path and not os.path.isabs(path):
                    restored[key] = os.path.abspath(
                        os.path.join(self.project_folder, path)
                    )
            ready[kind] = restored
        self._categorical_ready_configs = ready

        lai = state.get("lai_netcdf", {})
        self._lai_input_config = dict(lai) if isinstance(lai, dict) else {}
        path = self._lai_input_config.get("input_path")
        if path and not os.path.isabs(path):
            self._lai_input_config["input_path"] = os.path.abspath(
                os.path.join(self.project_folder, path)
            )
        self.task_coordinator.set_max_threads(state.get("thread_count", 2))

        for kind, text in state.get("categorical_types", {}).items():
            if kind not in {"lc", "soil", "geology", "lai"}:
                continue
            combo = self.categorical_type_combo(kind)
            index = combo.findText(str(text))
            if index < 0:
                normalized = str(text).lower().replace("_", " ")
                for candidate in range(combo.count()):
                    if (
                            combo.itemText(candidate).lower().replace("_", " ")
                            == normalized):
                        index = candidate
                        break
            combo.setCurrentIndex(index)
            self._categorical_modes[kind] = str(text) if index >= 0 else ""

    def restore_advanced_inputs(self, state):
        """Restore advanced source records without reopening their dialogs."""
        inputs = state.get("advanced_inputs", {})
        self._advanced_inputs = dict(inputs) if isinstance(inputs, dict) else {}
        source = str(state.get("land_cover_ready_source", "") or "")
        if source and not os.path.isabs(source):
            source = os.path.abspath(os.path.join(self.project_folder, source))
        self._land_cover_ready_source = source

    def restore_mhm_version(self, version):
        """Restore the saved mHM version selection."""
        combo_box = self.comboBox_mHMversion
        if combo_box is None or not version:
            return
        index = combo_box.findText(str(version))
        if index >= 0:
            combo_box.setCurrentIndex(index)

    def restore_grid_resolution_preferences(self, grid_resolutions):
        """Restore preferred L1/L11 selections for the next combo population."""
        try:
            self._preferred_l1_resolution = (
                float(grid_resolutions.get("l1_resolution"))
                if grid_resolutions.get("l1_resolution") is not None else None
            )
        except (TypeError, ValueError):
            self._preferred_l1_resolution = None
        try:
            self._preferred_l11_resolution = (
                float(grid_resolutions.get("l11_resolution"))
                if grid_resolutions.get("l11_resolution") is not None else None
            )
        except (TypeError, ValueError):
            self._preferred_l11_resolution = None

    def restore_text_inputs(self, text_inputs):
        """Restore saved file/folder and numeric text fields."""
        widget_lookup = dict(self.input_text_widgets())
        for key, value in text_inputs.items():
            widget = widget_lookup.get(key)
            if widget is not None:
                widget.setText(value or "")

    def restore_input_crs(self, crs_authid):
        """Restore the saved processing CRS if it is valid."""
        if not crs_authid:
            return

        from qgis.core import QgsCoordinateReferenceSystem

        crs = QgsCoordinateReferenceSystem(crs_authid)
        if crs.isValid():
            self.mProjectionSelectionWidget_crs.setCrs(crs)

    def restore_input_layers(self, saved_layers):
        """Restore saved layer selections, loading sources into QGIS if needed."""
        widget_lookup = dict(self.input_layer_widgets())
        for key, layer_info in saved_layers.items():
            combo_box = widget_lookup.get(key)
            if combo_box is None or not layer_info:
                continue
            if (
                isinstance(combo_box, InputComboAdapter)
                and layer_info.get("origin") == "file"
            ):
                self._restore_file_input(combo_box, layer_info)
                continue

            layer = self.find_or_load_saved_layer(layer_info)
            if not layer:
                self.log_message(
                    f"WARNING: Could not restore saved layer '{key}': "
                    f"{layer_info.get('source', '')}"
                )
                continue

            try:
                combo_box.setLayer(layer)
            except Exception as e:
                self.log_message(f"WARNING: Could not set saved layer '{key}': {e}")

    def _restore_file_input(self, adapter, layer_info):
        source = os.path.abspath(layer_info.get("source", ""))
        if not os.path.isfile(source):
            self.log_message(f"WARNING: Saved input file is missing: {source}")
            return
        combo = adapter.combo_box
        for index in range(combo.count()):
            data = combo.itemData(index)
            if isinstance(data, dict) and data.get("path") == source:
                combo.setCurrentIndex(index)
                return
        label = source
        try:
            relative = os.path.relpath(source, self.project_folder)
            if relative != ".." and not relative.startswith(f"..{os.sep}"):
                label = relative.replace("\\", "/")
        except ValueError:
            pass
        combo.addItem(
            label,
            {
                "origin": "file",
                "kind": adapter.kind,
                "path": source,
                "label": label,
                "manual": bool(layer_info.get("manual")),
            },
        )
        combo.setCurrentIndex(combo.count() - 1)

    def find_or_load_saved_layer(self, layer_info):
        """Find an existing QGIS layer by source or load it from disk."""
        source = layer_info.get("source", "")
        if not source:
            return None

        existing_layer = self.find_project_layer_by_source(source)
        if existing_layer:
            return existing_layer

        name = layer_info.get("name") or os.path.basename(source.split("|")[0])
        layer_type = layer_info.get("type", "vector")
        if layer_type == "raster":
            layer = QgsRasterLayer(source, name)
        else:
            layer = QgsVectorLayer(source, name, "ogr")

        if not layer.isValid():
            return None

        QgsProject.instance().addMapLayer(layer)
        self.log_message(f"Restored input layer: {name} | {source}")
        return layer

    def find_project_layer_by_source(self, source):
        """Find a loaded layer with the same source path/URI."""
        target_source = self.normalized_layer_source(source)
        for layer in QgsProject.instance().mapLayers().values():
            if self.normalized_layer_source(layer.source()) == target_source:
                return layer
        return None

    def normalized_layer_source(self, source):
        """Normalize a QGIS layer source for comparison."""
        if not source:
            return ""

        base_source = source.split("|")[0]
        if os.path.exists(base_source):
            base_source = os.path.abspath(base_source)
            suffix = source[len(source.split("|")[0]) :]
            return (base_source + suffix).replace("\\", "/").lower()

        return source.replace("\\", "/").lower()

    def log_selected_input_paths(self, action_name):
        """Print selected layer/file inputs for easier debugging from QGIS."""
        self.log_message(f"\n--- Input selections before {action_name} ---")
        self.log_message(f"Project folder: {self.project_folder or '<not selected>'}")

        layer_inputs = [
            ("DEM", self.input_combo("dem")),
            ("Pour points", self.input_combo("pour_points")),
            ("Land cover", self.input_combo("land_cover")),
            ("Soil", self.input_combo("soil")),
            ("Geology", self.input_combo("geology")),
        ]

        layer_inputs.append(("LAI class", self.input_combo("lai")))

        for label, combo_box in layer_inputs:
            layer = combo_box.currentLayer()
            if layer:
                self.log_message(f"{label}: {layer.name()} | {layer.source()}")
            else:
                self.log_message(f"{label}: <not selected>")

        for kind, label in (
            ("lc", "Land cover"),
            ("soil", "Soil"),
            ("geology", "Geology"),
        ):
            mode = self.categorical_input_mode(kind) or "<not selected>"
            self.log_message(f"{label} input type: {mode}")
            config = self.categorical_lookup_config(kind)
            if config and mode.lower() == "lookup table":
                self.log_message(
                    f"{label} lookup: {config['lookup_table']} | "
                    f"{config['mapping_field']} -> {config['class_field']}"
                )

        self.log_message(
            f"LAI input type: {self.comboBox_lai_inputType.currentText() or '<not selected>'}"
        )
        if hasattr(self, "lineEdit_lai_file"):
            self.log_message(
                f"LAI file: {self.lineEdit_lai_file.text() or '<not selected>'}"
            )
        for kind, _, _ in self.meteo_input_widgets():
            self.log_message(
                f"{kind.title()} folder: "
                f"{self.selected_meteo_folder(kind) or '<not selected>'}"
            )
            self.log_message(
                f"{kind.title()} data source: "
                f"{self.selected_meteo_source(kind) or '<not selected>'}"
            )

    # --- Project Management Methods ---

    def select_project_folder(self):
        """Opens a dialog to select the project working directory."""
        folder = QFileDialog.getExistingDirectory(self, "Select Project Folder")
        if folder:
            self.project_folder = folder
            self.lineEdit_ProjectFolder.setText(self.project_folder)
            self.log_message(f"Project folder set to: {self.project_folder}")

            created = ensure_project_structure(
                self.project_folder, self.configuration_processor.selected_version()
            )
            self.geometry_folder = geometry_folder(self.project_folder)
            if created:
                self.log_message(
                    f"Project structure prepared with {len(created)} folder(s)."
                )

            self.refresh_input_sources()
            self.refresh_meteo_folder_sources()
            self.refresh_meteo_display()
            self.load_input_state()
            self.morphology_processor.update_gauged_outlet_count()

            # Load project state in morphology processor
            self.morphology_processor.load_project_state()
            self.refresh_morphology_workflow_button_states()
            self.meteorology_processor.load_project_state()
            self.refresh_grid_resolution_controls()

    def on_tab_changed(self, index):
        """Switches the stacked widget page when the tab is changed."""
        page = self.page_for_tab_index(index)
        if page is not None:
            self.stackedWidget.setCurrentWidget(page)
        else:
            self.stackedWidget.setCurrentIndex(index)
        self.log_message(f"Switched to '{self.tabWidget.tabText(index)}' tab.")

    def page_for_tab_index(self, index):
        """Return the stacked page associated with a tab widget index."""
        tab = self.tabWidget.widget(index)
        page_pairs = (
            (self.tab_geometry, self.page_geometry),
            (self.tab_meteo, self.page_meteo),
            (self.tab_hydro, self.page_hydro),
            (
                self.tab_configuration,
                self.page_configuration,
            ),
            (
                self.tab_execution,
                self.page_execution,
            ),
            (
                self.tab_calibration,
                self.page_calibration,
            ),
            (self.tab_outputs, self.page_outputs),
        )
        for tab_widget, page_widget in page_pairs:
            if tab_widget is tab:
                return page_widget
        return None

    # --- UI Helper Methods ---

    def browse_lai_file(self):
        """Browse for Leaf Area Index file"""
        if not hasattr(self, "lineEdit_lai_file"):
            self.log_message("LAI input is now selected from the layer dropdown.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Leaf Area Index File",
            "",
            "All Files (*);;GeoTIFF (*.tif *.tiff);;NetCDF (*.nc);;HDF5 (*.h5 *.hdf5)",
        )
        if file_path:
            self.lineEdit_lai_file.setText(file_path)
            self.log_message(f"LAI file selected: {file_path}")
            self.save_input_state()
            self.invalidate_meteo_morph_setup()

    def get_lai_time_range(self):
        """Get the selected LAI time/date for extraction"""
        if hasattr(self, "dateEdit") and self.dateEdit:
            return self.dateEdit.dateTime()
        return None

    def get_crs(self):
        """Get the selected CRS"""
        return self.mProjectionSelectionWidget_crs.crs()

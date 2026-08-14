# -*- coding: utf-8 -*-
"""Dialog and UI wiring for the pymhm QGIS plugin."""

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
    from qgis.PyQt.QtWidgets import QDialog, QFileDialog, QMessageBox
except ImportError:
    from .standalone_qgis import install

    install(force=True)
    from qgis.core import (QgsApplication, QgsMapLayer, QgsProject,
                           QgsRasterLayer, QgsVectorLayer)
    from qgis.PyQt import QtCore
    from qgis.PyQt.QtWidgets import QDialog, QFileDialog, QMessageBox

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
                              is_geographic_unit, load_meteo_grid_metadata,
                              possible_resolutions, raster_resolution_info,
                              read_header_file)
from .input_selection import (INPUT_EXTENSIONS, InputComboAdapter,
                              LookupConfigDialog, loaded_qgis_items,
                              scan_project_folders, scan_project_inputs)
from .Meteorology import MeteorologyProcessor
from .Meteorology.forcing import (MeteoFolderSpec, TargetGrid,
                                  inspect_meteo_folder,
                                  inspect_meteo_inputs, resolution_in_crs)
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
from .morphology_display import (
    DISPLAY_KEYS,
    land_cover_periods,
    resolve_display_output,
)
from .qgis_compat import map_layer_filters
from .terminal_dialog import ProjectTerminalDialog
from .pyui.ui_pymhm_dialog_base import Ui_pymhmDialog
# Import utility mixin and processors
from .utils import DialogUtils


class MorphologyWorkflowWorker(QtCore.QObject):
    """Run a morphology workflow away from the dialog event loop."""

    log_message = QtCore.pyqtSignal(str)
    finished = QtCore.pyqtSignal(str, bool, str)

    def __init__(
            self,
            processor,
            workflow_key,
            workflow_label,
            method_name,
            project_folder,
            dem_layer,
            pour_points_layer):
        super().__init__()
        self.processor = processor
        self.workflow_key = workflow_key
        self.workflow_label = workflow_label
        self.method_name = method_name
        self.project_folder = project_folder
        self.dem_layer = dem_layer
        self.pour_points_layer = pour_points_layer
        self._original_log_message = None
        self._original_run_processing_algorithm = None
        self._original_check_prerequisites = None

    @QtCore.pyqtSlot()
    def run(self):
        """Execute the workflow and emit the result."""
        self._install_worker_hooks()
        try:
            workflow_method = getattr(self.processor, self.method_name)
            ok = bool(workflow_method(show_error_dialog=False))
            self.finished.emit(self.workflow_key, ok, "")
        except Exception as exc:
            details = traceback.format_exc()
            self.log_message.emit(
                f"\nERROR: {self.workflow_label} worker failed with exception: {exc}"
            )
            self.log_message.emit(f"Traceback: {details}")
            try:
                self.processor.mark_workflow_status(
                    self.workflow_key,
                    "failed",
                    f"{self.workflow_label} worker failed: {exc}",
                )
            except Exception:
                pass
            self.finished.emit(self.workflow_key, False, str(exc))
        finally:
            self._restore_worker_hooks()

    def _install_worker_hooks(self):
        self._original_log_message = self.processor.log_message
        self._original_run_processing_algorithm = (
            self.processor.run_processing_algorithm
        )
        self._original_check_prerequisites = self.processor.check_prerequisites
        self.processor.log_message = self.log_message.emit
        self.processor.run_processing_algorithm = self._run_processing_algorithm
        self.processor.check_prerequisites = self._check_prerequisites

    def _restore_worker_hooks(self):
        if self._original_log_message is not None:
            self.processor.log_message = self._original_log_message
        if self._original_run_processing_algorithm is not None:
            self.processor.run_processing_algorithm = (
                self._original_run_processing_algorithm
            )
        if self._original_check_prerequisites is not None:
            self.processor.check_prerequisites = self._original_check_prerequisites

    def _check_prerequisites(self, needs_pour_points=False):
        if not self.project_folder:
            self.log_message.emit(
                "ERROR: Please select a project folder before proceeding."
            )
            return False
        if not self.dem_layer:
            self.log_message.emit("ERROR: Please select a DEM Raster Layer.")
            return False
        if needs_pour_points and not self.pour_points_layer:
            self.log_message.emit("ERROR: Please select a Pour Points Layer.")
            return False
        return True

    def _run_processing_algorithm(self, name, params):
        import processing

        self.log_message.emit(f"Running algorithm: {name}...")
        try:
            result = processing.run(name, params)
            self.log_message.emit(f"Algorithm '{name}' finished successfully.")
            self.processor.record_processing_outputs(name, params, result)
            return result
        except Exception as exc:
            self.log_message.emit(
                f"ERROR: Algorithm '{name}' failed. Details: {exc}"
            )
            return None


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


class pymhmDialog(QDialog, Ui_pymhmDialog, DialogUtils):
    def __init__(self, parent=None):
        """Constructor."""

        super(pymhmDialog, self).__init__(parent)
        self.setupUi(self)
        self._categorical_lookup_configs = {}
        self._categorical_modes = {}
        self._advanced_inputs = {}
        self._land_cover_ready_source = ""
        self._domain_definition_mode = ""
        self._input_adapters = {}
        self.configure_widget_aliases()
        self.configure_morphology_display()

        # --- Filter map layer combo boxes to show only relevant layer types ---
        self.mMapLayerComboBox_dem.setFilters(map_layer_filters("RasterLayer"))
        self.mMapLayerComboBox_pour_points.setFilters(map_layer_filters("VectorLayer"))

        # Set filters for new layer combo boxes (both vector and raster allowed)
        self.mMapLayerComboBox_soil.setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        if hasattr(self, "mMapLayerComboBox_land_cover"):
            self.mMapLayerComboBox_land_cover.setFilters(
                map_layer_filters("RasterLayer", "VectorLayer")
            )
        self.mMapLayerComboBox_geology.setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            self.mMapLayerComboBox_LAI_Class.setFilters(
                map_layer_filters("RasterLayer", "VectorLayer")
            )
        self.configure_input_layer_combo_boxes()

        # --- Instance attributes for managing file paths ---
        self.project_folder = None
        self.geometry_folder = None  # Subfolder for geometry outputs
        self.input_state_filename = "pymhm_input_state.json"
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

        # --- Connect signals and slots ---
        self.configure_page_aliases()
        self.connect_signals()
        self.refresh_input_sources()
        self.refresh_meteo_folder_sources()
        self.refresh_grid_resolution_controls()

    def configure_widget_aliases(self):
        """Keep renamed widgets compatible with older dialog code."""

        input_widgets = {
            "dem": ("comboBox_demInput", "mMapLayerComboBox_dem"),
            "pour_points": (
                "comboBox_pourPointInput",
                "mMapLayerComboBox_pour_points",
            ),
            "land_cover": (
                "comboBox_landUseInput",
                "mMapLayerComboBox_land_cover",
            ),
            "soil": ("comboBox_soilInput", "mMapLayerComboBox_soil"),
            "geology": (
                "comboBox_geologyInput",
                "mMapLayerComboBox_geology",
            ),
            "lai": ("comboBox_laiInput", "mMapLayerComboBox_LAI_Class"),
        }
        for kind, (combo_name, legacy_name) in input_widgets.items():
            combo = getattr(self, combo_name, None)
            if combo is None:
                continue
            adapter = InputComboAdapter(combo, kind, self)
            self._input_adapters[kind] = adapter
            setattr(self, legacy_name, adapter)

        if hasattr(self, "comboBox_lai_inputType"):
            self.comboBox_laiInputType = self.comboBox_lai_inputType
        if (
            hasattr(self, "comboBox_landUseInputType")
            and not hasattr(self, "comboBox_landCover_inputType")
        ):
            self.comboBox_landCover_inputType = self.comboBox_landUseInputType

        if (
            hasattr(self, "pushButton_executeAllMorphology")
            and not hasattr(self, "pushButton_executeAll")
        ):
            self.pushButton_executeAll = self.pushButton_executeAllMorphology
        if (
            hasattr(self, "pushButton_execute_mHM")
            and not hasattr(self, "pushButton_RUN")
        ):
            self.pushButton_RUN = self.pushButton_execute_mHM

        for current_name, folder_name in (
            ("pushButton_browsePrecipitationFile", "pushButton_browsePrecipitationFolder"),
            ("pushButton_browseTemperatureFile", "pushButton_browseTemperatureFolder"),
            ("pushButton_browsePetFile", "pushButton_browsePetFolder"),
        ):
            current = getattr(self, current_name, None)
            if current is not None and not hasattr(self, folder_name):
                setattr(self, folder_name, current)

    def configure_page_aliases(self):
        """Keep renamed stacked pages compatible with older plugin code."""

        if hasattr(self, "page_configuration") and not hasattr(self, "page_validation"):
            self.page_validation = self.page_configuration
        if hasattr(self, "page_execution") and not hasattr(self, "page_datasets"):
            self.page_datasets = self.page_execution
        if hasattr(self, "pushButton_Terminal") and not hasattr(self, "pushButton_terminal"):
            self.pushButton_terminal = self.pushButton_Terminal

    def configure_morphology_display(self):
        """Attach stable keys to the morphology display choices."""
        combo = getattr(self, "comboBox_morphVariableToDisplay", None)
        if combo is not None:
            for index, key in enumerate(DISPLAY_KEYS[: combo.count()]):
                combo.setItemData(index, key)
        editor = getattr(self, "dateTimeEdit_forHistoricalInputs", None)
        if editor is not None:
            editor.setEnabled(False)

    def configure_input_layer_combo_boxes(self):
        """Allow input layer boxes to start empty so layers are chosen deliberately."""

        layer_combo_boxes = [
            combo
            for combo in (
                getattr(self, "mMapLayerComboBox_dem", None),
                getattr(self, "mMapLayerComboBox_pour_points", None),
                getattr(self, "mMapLayerComboBox_soil", None),
                getattr(self, "mMapLayerComboBox_land_cover", None),
                getattr(self, "mMapLayerComboBox_geology", None),
            )
            if combo is not None
        ]

        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            layer_combo_boxes.append(self.mMapLayerComboBox_LAI_Class)

        for combo_box in layer_combo_boxes:
            if hasattr(combo_box, "setAllowEmptyLayer"):
                try:
                    combo_box.setAllowEmptyLayer(True)
                except TypeError:
                    combo_box.setAllowEmptyLayer(True, "")

            if hasattr(combo_box, "setLayer"):
                try:
                    combo_box.setLayer(None)
                except TypeError:
                    combo_box.setCurrentIndex(-1)
            else:
                combo_box.setCurrentIndex(-1)

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
            layer = self.mMapLayerComboBox_pour_points.currentLayer()

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
        layer = self.mMapLayerComboBox_pour_points.currentLayer()
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
            metadata = inspect_meteo_folder(
                folder,
                kind,
                source,
                crs_fallback=(
                    self.selected_meteo_crs()
                    if source == "mhm_ready" else None
                ),
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

    def connect_input_source_signals(self):
        """Connect folder searching, browsing, and categorical type controls."""
        self.checkBox_enableFolderSearch.toggled.connect(
            lambda checked=False: (
                self.refresh_input_sources(),
                self.save_input_state(),
            )
        )
        for kind, button_name in (
            ("dem", "pushButton_browseDEMInput"),
            ("pour_points", "pushButton_browsePourPointInput"),
            ("land_cover", "pushButton_browseLandUseInput"),
            ("soil", "pushButton_soilInput"),
            ("geology", "pushButton_geologyInput"),
            ("lai", "pushButton_laiInput"),
        ):
            button = getattr(self, button_name, None)
            if button is not None:
                button.setToolTip(
                    f"Browse for {kind.replace('_', ' ')} input file"
                )
                button.clicked.connect(
                    lambda checked=False, input_kind=kind: self.browse_input_file(
                        input_kind
                    )
                )

        for kind, combo_name in (
            ("lc", "comboBox_landCover_inputType"),
            ("soil", "comboBox_soil_inputType"),
            ("geology", "comboBox_geology_inputType"),
        ):
            combo = getattr(self, combo_name, None)
            if combo is not None:
                combo.activated.connect(
                    lambda _index, input_kind=kind, widget=combo: (
                        self.handle_categorical_type(
                            input_kind, widget.currentText()
                        )
                    )
                )

        project = QgsProject.instance()
        for signal_name in ("layersAdded", "layersRemoved"):
            signal = getattr(project, signal_name, None)
            if signal is not None:
                try:
                    signal.connect(lambda *_args: self.refresh_input_sources())
                except Exception:
                    pass

    def connect_signals(self):
        """Connect all UI element signals to appropriate slots."""

        # Project management
        self.pushButton_BrowseProjectFolder.clicked.connect(self.select_project_folder)
        self.tabWidget.currentChanged.connect(self.on_tab_changed)
        self.connect_input_source_signals()

        # Morphology/Geometry processing - delegate to processor
        self.connect_optional_processor_button(
            "pushButton_convertDEMtoASC",
            "Convert DEM to ASC",
            self.morphology_processor.convert_dem_to_asc,
        )
        self.connect_optional_processor_button(
            "pushButton_fillDem", "Fill DEM", self.morphology_processor.fill_dem
        )
        self.connect_optional_processor_button(
            "pushButton_createNetwork",
            "Create Channel Network",
            self.morphology_processor.process_channel_network,
        )
        self.connect_optional_processor_button(
            "pushButton_snapPoints",
            "Snap Pour Points",
            self.morphology_processor.snap_points,
        )
        if hasattr(self, "pushButton_delineate"):
            self.pushButton_delineate.clicked.connect(self.open_domain_delineator)
        if hasattr(self, "pushButton_elevation_bands"):
            self.connect_processor_button(
                self.pushButton_elevation_bands,
                "Elevation Bands",
                self.morphology_processor.process_elevation_bands,
            )
        if hasattr(self, "pushButton_bandDetails"):
            self.connect_processor_button(
                self.pushButton_bandDetails,
                "Elevation Band Land Cover Details",
                self.morphology_processor.process_band_details,
            )

        # Hydrological processing - delegate to processor
        self.connect_optional_processor_button(
            "pushButton_aspect", "Aspect", self.morphology_processor.process_aspect
        )
        self.connect_optional_processor_button(
            "pushButton_slope", "Slope", self.morphology_processor.process_slope
        )
        self.connect_optional_processor_button(
            "pushButton_flowAccumulation",
            "Flow Accumulation",
            self.morphology_processor.process_flow_accumulation,
        )
        self.connect_optional_processor_button(
            "pushButton_flowDirection",
            "Flow Direction",
            self.morphology_processor.process_flow_direction,
        )
        self.connect_optional_processor_button(
            "pushButton_gaugePosition",
            "Gauge Position",
            self.morphology_processor.process_gauge_position,
        )
        if hasattr(self, "pushButton_assignDischargeTables"):
            self.connect_processor_button(
                self.pushButton_assignDischargeTables,
                "Assign Discharge Tables",
                self.morphology_processor.assign_discharge_tables,
            )
        try:
            self.mMapLayerComboBox_pour_points.layerChanged.connect(
                self.populate_pour_point_outlet_fields
            )
            self.mMapLayerComboBox_pour_points.layerChanged.connect(
                lambda layer=None: self.morphology_processor.update_gauged_outlet_count(
                    layer
                )
            )
        except Exception:
            pass
        self.morphology_processor.update_gauged_outlet_count()

        self.comboBox_pourPointOutletID.currentIndexChanged.connect(
            self.handle_model_input_changed
        )
        self.checkBox_DEMdomain.toggled.connect(
            self.handle_model_input_changed
        )
        if hasattr(self, "comboBox_domainDefinitionType"):
            self.comboBox_domainDefinitionType.activated.connect(
                self.handle_domain_definition_type
            )

        # Layer processing - delegate to processor
        self.connect_optional_processor_button(
            "pushButton_landUse",
            "Land Use",
            self.morphology_processor.process_land_use,
        )
        self.connect_optional_processor_button(
            "pushButton_soil", "Soil", self.morphology_processor.process_soil
        )
        self.connect_optional_processor_button(
            "pushButton_hydrogeology",
            "Hydrogeology",
            self.morphology_processor.process_geology,
        )
        if hasattr(self, "pushButton_LAI"):
            self.connect_processor_button(
                self.pushButton_LAI,
                "LAI",
                self.morphology_processor.process_lai,
            )
        if hasattr(self, "pushButton_morphLayerDisplay"):
            self.pushButton_morphLayerDisplay.clicked.connect(
                self.display_selected_morphology_layer
            )
        if hasattr(self, "comboBox_morphVariableToDisplay"):
            self.comboBox_morphVariableToDisplay.currentIndexChanged.connect(
                self.update_morphology_date_control
            )

        # Reset geometry
        self.pushButton_resetGeometry.clicked.connect(self.reset_geometry_processing)

        # Batch morphology workflows
        execute_all_button = self.morphology_workflow_button("execute_all")
        if execute_all_button is not None:
            execute_all_button.clicked.connect(self.start_execute_all_processing)
        self.pushButton_executeMeteoMorphSetup.clicked.connect(
            self.start_meteo_morph_setup_processing
        )

        # LAI file browser
        if hasattr(self, "pushButton_browse_lai"):
            self.pushButton_browse_lai.clicked.connect(self.browse_lai_file)

        self.connect_grid_resolution_signals()

        for kind, button in (
            ("precipitation", self.pushButton_browsePrecipitationFolder),
            ("temperature", self.pushButton_browseTemperatureFolder),
            ("pet", self.pushButton_browsePetFolder),
        ):
            button.clicked.connect(
                lambda checked=False, meteo_kind=kind:
                self.browse_meteo_input_folder(meteo_kind)
            )
        for kind, folder_combo, source_combo in self.meteo_input_widgets():
            folder_combo.currentIndexChanged.connect(
                lambda index=None, meteo_kind=kind:
                self.handle_meteo_input_changed(meteo_kind)
            )
            source_combo.currentIndexChanged.connect(
                lambda index=None, meteo_kind=kind:
                self.handle_meteo_input_changed(meteo_kind)
            )

        # Configuration/execution processing - delegate to processor
        if hasattr(self, "pushButton_createLatLon"):
            self.connect_processor_button(
                self.pushButton_createLatLon,
                "Create LatLon",
                self.morphology_processor.process_lat_lon,
            )
            self.update_latlon_button_state()
        if hasattr(self, "pushButton_edit_nmls"):
            self.pushButton_edit_nmls.setToolTip("Edit mHM namelists")
            self.connect_processor_button(
                self.pushButton_edit_nmls,
                "Edit Namelists",
                self.configuration_processor.edit_namelists,
            )
        run_mhm_button = getattr(self, "pushButton_execute_mHM", None)
        if run_mhm_button is None:
            run_mhm_button = getattr(self, "pushButton_RUN", None)
        if run_mhm_button is not None:
            run_mhm_button.setToolTip("Run mHM")
            self.connect_processor_button(
                run_mhm_button,
                "Run mHM",
                self.configuration_processor.run_mhm,
            )
        connected_terminal_buttons = set()
        for terminal_button_name in ("pushButton_terminal", "pushButton_Terminal"):
            terminal_button = getattr(self, terminal_button_name, None)
            if (
                terminal_button is not None
                and id(terminal_button) not in connected_terminal_buttons
            ):
                terminal_button.clicked.connect(self.open_project_terminal)
                connected_terminal_buttons.add(id(terminal_button))
        if hasattr(self, "comboBox_mHMversion"):
            self.comboBox_mHMversion.currentIndexChanged.connect(
                lambda index=None: self.configuration_processor.handle_version_changed()
            )
        self.connect_input_state_signals()
        self.connect_input_state_teardown_guards()

        # Initialize CRS widget with project CRS
        project_crs = QgsProject.instance().crs()
        if project_crs.isValid():
            self.mProjectionSelectionWidget_crs.setCrs(project_crs)

    def connect_optional_processor_button(self, name, label, callback):
        """Connect a processing control when it exists in the active UI."""
        button = getattr(self, name, None)
        if button is not None:
            self.connect_processor_button(button, label, callback)

    def connect_grid_resolution_signals(self):
        """Refresh derived grid controls when DEM/L1/L11 selections change."""
        try:
            self.mMapLayerComboBox_dem.layerChanged.connect(
                lambda layer=None: self.update_l0_resolution_from_dem(layer)
            )
            self.mMapLayerComboBox_dem.layerChanged.connect(
                lambda layer=None:
                self.inspect_meteo_selection("precipitation", show_errors=False)
            )
        except Exception:
            pass

        self.spinBox_L2ResolutionMultiplier.valueChanged.connect(
            self.handle_l2_multiplier_changed
        )
        try:
            self.mProjectionSelectionWidget_crs.crsChanged.connect(
                lambda crs=None:
                self.inspect_meteo_selection("precipitation", show_errors=False)
            )
        except Exception:
            pass

        if hasattr(self, "comboBox_L1"):
            try:
                self.comboBox_L1.currentIndexChanged.connect(
                    lambda index=None: self.handle_l1_resolution_changed()
                )
            except Exception:
                pass

        if hasattr(self, "comboBox_L11"):
            try:
                self.comboBox_L11.currentIndexChanged.connect(
                    lambda index=None: self.update_l11_resolution_label()
                )
            except Exception:
                pass

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
        if not self.mMapLayerComboBox_dem.currentLayer():
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
                header = dict(header)
                header["cellsize"] = ceil_cellsize(header["cellsize"], unit)
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
            getattr(self, "comboBox_L1", None),
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
        if hasattr(self, "label_L1ResolutionUnit"):
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
            getattr(self, "comboBox_L11", None),
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
        if hasattr(self, "label_L11ResolutionUnit"):
            self.label_L11ResolutionUnit.setText(unit if l11_values else "")
        self.update_l11_resolution_label()
        self.update_latlon_button_state()

    def disable_l1_l11_resolution_options(self):
        """Disable L1/L11 controls until meteo L2 grid metadata exists."""
        self.disable_resolution_combo(getattr(self, "comboBox_L1", None))
        self.disable_l11_resolution_options()
        self._set_resolution_labels("L1", "", "")

    def disable_l11_resolution_options(self):
        """Disable L11 controls and clear its labels."""
        self.disable_resolution_combo(getattr(self, "comboBox_L11", None))
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
        label = getattr(self, "label_L1Resolution", None)
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
        label = getattr(self, "label_L11Resolution", None)
        if label is None:
            return
        if value and l1_resolution:
            label.setText(f"{value / l1_resolution:g} x L1")
        else:
            label.setText("")
        self.update_latlon_button_state()

    def current_l0_resolution(self):
        """Return current L0 resolution."""
        if self._grid_l0_info:
            return float(self._grid_l0_info["resolution"])
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
        return self._current_combo_resolution(getattr(self, "comboBox_L1", None))

    def current_l11_resolution(self):
        """Return selected L11 resolution."""
        return self._current_combo_resolution(getattr(self, "comboBox_L11", None))

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
        l2_header["cellsize"] = ceil_cellsize(l2_header["cellsize"], unit)
        l2_header["unit"] = unit
        return {
            "L0": header_for_existing_bounds(l2_header, l0_resolution, unit),
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
        button = getattr(self, "pushButton_createLatLon", None)
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

    def connect_input_state_signals(self):
        """Save input selections when editable inputs change."""

        for _, combo_box in self.input_layer_widgets():
            try:
                combo_box.layerChanged.connect(
                    self.handle_model_input_changed
                )
            except Exception:
                pass

        for _, line_edit in self.input_text_widgets():
            try:
                line_edit.editingFinished.connect(self.handle_model_input_changed)
            except Exception:
                try:
                    line_edit.textChanged.connect(
                        self.handle_model_input_changed
                    )
                except Exception:
                    pass

        if hasattr(self, "comboBox_laiInputType"):
            try:
                self.comboBox_laiInputType.currentIndexChanged.connect(
                    self.handle_model_input_changed
                )
            except Exception:
                pass

        for combo_box in (
                getattr(self, "comboBox_L1", None),
                getattr(self, "comboBox_L11", None)):
            if combo_box is None:
                continue
            try:
                combo_box.currentIndexChanged.connect(
                    self.handle_model_input_changed
                )
            except Exception:
                pass

        try:
            self.mProjectionSelectionWidget_crs.crsChanged.connect(
                self.handle_model_input_changed
            )
        except Exception:
            pass

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

    def connect_processor_button(self, button, action_name, callback):
        """Connect a button to a processor callback with input path logging."""
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
        if getattr(sys.modules.get("qgis"), "_pymhm_standalone_qgis", False):
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

            layer = self.mMapLayerComboBox_pour_points.currentLayer()
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

    def update_morphology_date_control(self, index=None):
        """Enable the date selector only for historical land cover."""
        editor = getattr(self, "dateTimeEdit_forHistoricalInputs", None)
        combo = getattr(self, "comboBox_morphVariableToDisplay", None)
        if editor is None or combo is None:
            return
        key = combo.currentData()
        periods = land_cover_periods(self.project_folder) if self.project_folder else []
        enabled = key == "land_cover" and len(periods) > 1
        editor.setEnabled(enabled)
        if not enabled:
            return
        try:
            first = min(int(period["start_year"]) for period in periods)
            last = max(int(period["end_year"]) for period in periods)
        except (KeyError, TypeError, ValueError):
            editor.setEnabled(False)
            return
        editor.setMinimumDate(QtCore.QDate(first, 1, 1))
        editor.setMaximumDate(QtCore.QDate(last, 12, 31))
        editor.setDate(QtCore.QDate(first, 1, 1))

    def display_selected_morphology_layer(self, checked=False):
        """Load the prepared output selected by the display combo box."""
        if not self.project_folder:
            QMessageBox.warning(self, "Display Layer", "Select a project folder first.")
            return
        combo = self.comboBox_morphVariableToDisplay
        key = combo.currentData()
        if not key:
            QMessageBox.warning(self, "Display Layer", "Select a morphology layer.")
            return
        year = None
        if self.dateTimeEdit_forHistoricalInputs.isEnabled():
            year = self.dateTimeEdit_forHistoricalInputs.date().year()
        output = resolve_display_output(self.project_folder, key, year=year)
        if output is None:
            QMessageBox.warning(
                self,
                "Display Layer",
                "The selected morphology output has not been prepared yet.",
            )
            return
        source = str(output.path)
        if output.variable and output.path.suffix.lower() == ".nc":
            source = f'NETCDF:"{output.path}":{output.variable}'
        layer = self.load_layer(source, output.name, output.is_raster)
        if layer is not None and output.band and output.band > 1:
            try:
                from qgis.core import QgsSingleBandGrayRenderer

                layer.setRenderer(
                    QgsSingleBandGrayRenderer(layer.dataProvider(), output.band)
                )
                layer.triggerRepaint()
            except Exception as error:
                self.log_message(
                    f"WARNING: Could not select historical band {output.band}: {error}"
                )

    def open_domain_delineator(self, checked=False):
        """Open the per-outlet domain delineation dialog in QGIS."""
        if getattr(sys.modules.get("qgis"), "_pymhm_standalone_qgis", False):
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

        try:
            from .domain_delineator_dialog import DomainDelineatorDialog

            self.morphology_processor.load_project_state()
            dialog = DomainDelineatorDialog(
                self,
                self.morphology_processor,
                self.mMapLayerComboBox_pour_points.currentLayer(),
                self.selected_outlet_id_field(),
                outlet_ids,
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

    def morphology_workflow_specs(self):
        """Return metadata for threaded morphology workflows."""
        return {
            "execute_all": {
                "button": "pushButton_executeAllMorphology",
                "fallback_button": "pushButton_executeAll",
                "label": "Execute All Processing",
                "action_name": "Execute All Processing",
                "method": "execute_all_processing",
                "thread_name": "PymHMExecuteAllThread",
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
                "thread_name": "PymHMMeteoMorphSetupThread",
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
        for workflow_key, thread in self._morphology_workflow_threads.items():
            if thread is not None and thread.isRunning():
                return workflow_key
        return None

    def start_execute_all_processing(self):
        """Run execute-all processing in a background worker thread."""
        self.start_morphology_workflow("execute_all")

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
        spec = self.morphology_workflow_specs().get(workflow_key)
        if spec is None:
            self.log_message(f"ERROR: Unknown morphology workflow: {workflow_key}")
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
        worker_arguments = (
            self.morphology_processor,
            workflow_key,
            spec["label"],
            spec["method"],
            self.project_folder,
            self.mMapLayerComboBox_dem.currentLayer(),
            self.mMapLayerComboBox_pour_points.currentLayer(),
        )
        if workflow_key == "meteo_morph_setup":
            worker = MeteorologyWorkflowWorker(
                workflow_key,
                spec["label"],
                meteorology_processor=self.meteorology_processor,
                meteorology_run=self._pending_meteo_run,
            )
        else:
            worker = MorphologyWorkflowWorker(*worker_arguments)
        worker.moveToThread(thread)

        worker.log_message.connect(self.log_message)
        if workflow_key == "meteo_morph_setup":
            worker.finished.connect(self.finish_meteo_workflow_stage)
        else:
            worker.finished.connect(self.finish_morphology_workflow)
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
        if running_key is not None:
            spec = self.morphology_workflow_specs().get(running_key, {})
            QMessageBox.warning(
                self,
                spec.get("label", "Morphology Processing"),
                (
                    f"{spec.get('label', 'Morphology processing')} is still running. "
                    "Wait for it to finish before closing PymHM."
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
        }[kind]
        combo = getattr(self, name, None)
        if combo is None and kind == "lc":
            combo = getattr(self, "comboBox_landCover_inputType")
        return combo

    def categorical_input_mode(self, kind):
        """Return the selected categorical input mode."""
        return self.categorical_type_combo(kind).currentText().strip()

    def categorical_lookup_config(self, kind):
        """Return the accepted lookup-table selection for one data type."""
        config = self._categorical_lookup_configs.get(kind)
        return dict(config) if config else None

    def handle_categorical_type(self, kind, text):
        """Open the lookup dialog immediately when lookup mode is selected."""
        text = str(text or "").strip()
        previous = self._categorical_modes.get(kind, "")
        if self._loading_input_state:
            self._categorical_modes[kind] = text
            return
        normalized = text.lower()
        if kind == "lc":
            self.handle_land_use_input_type(text, previous)
            return
        if kind == "soil" and "multi-horizon" in normalized:
            self.handle_multi_horizon_soil_input(text, previous)
            return
        if "lookup table" not in normalized:
            if kind == "soil":
                self._advanced_inputs.pop("soil", None)
                self._save_standard_soil_nml_input("mhm_ready")
            self._categorical_modes[kind] = text
            self.save_input_state()
            self.invalidate_meteo_morph_setup()
            return
        if not self.project_folder:
            QMessageBox.warning(
                self,
                "Project Folder Required",
                "Select a project folder before configuring a lookup table.",
            )
            self._restore_categorical_mode(kind, previous)
            return

        dialog = LookupConfigDialog(
            self.project_folder,
            self,
            initial=self._categorical_lookup_configs.get(kind),
        )
        execute = getattr(dialog, "exec", None) or dialog.exec_
        if execute() != QDialog.Accepted or dialog.selected_config() is None:
            self._restore_categorical_mode(kind, previous)
            return

        config = dialog.selected_config()
        self._categorical_lookup_configs[kind] = {
            "lookup_table": config.lookup_table,
            "mapping_field": config.mapping_field,
            "class_field": config.class_field,
        }
        self._categorical_modes[kind] = text
        if kind == "soil":
            self._advanced_inputs.pop("soil", None)
            self._save_standard_soil_nml_input("single_categorical")
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
        if str(text).strip().lower() == "mhm ready":
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Select mHM-ready Land Cover",
                self.project_folder,
                "Land-cover rasters (*.asc *.nc *.tif);;All files (*)",
            )
            if not path:
                self._restore_categorical_mode("lc", previous)
                return
            try:
                from .advanced_input_processing import configure_ready_land_cover

                output = configure_ready_land_cover(
                    self.project_folder,
                    path,
                    self.comboBox_mHMversion.currentText(),
                )
                self._land_cover_ready_source = str(Path(path).resolve())
                self._advanced_inputs.pop("land_cover", None)
                self.morphology_processor.mark_output_prepared(
                    str(output), name=output.name, loaded=False
                )
                self.log_message(f"mHM-ready land cover saved: {output}")
            except Exception as error:
                QMessageBox.critical(self, "Land-use Input", str(error))
                self._restore_categorical_mode("lc", previous)
                return
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
        return kind == "lc" or (kind == "soil" and "multi-horizon" in text)

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
            else "data/static/morph/soil_class.asc"
        )
        definition = (
            relative_workspace_path(self.project_folder, classdefinition_path)
            if classdefinition_path
            else "data/static/morph/soil_classdefinition.txt"
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
        combo_box = getattr(self, "comboBox_laiInputType", None)
        if combo_box is None or not lai_input_type:
            return
        index = combo_box.findText(lai_input_type)
        if index >= 0:
            combo_box.setCurrentIndex(index)

    def input_layer_widgets(self):
        """Return persistent layer input widgets and their state keys."""
        widgets = []
        for key, name in (
            ("dem", "mMapLayerComboBox_dem"),
            ("pour_points", "mMapLayerComboBox_pour_points"),
            ("soil", "mMapLayerComboBox_soil"),
            ("land_cover", "mMapLayerComboBox_land_cover"),
            ("geology", "mMapLayerComboBox_geology"),
        ):
            widget = getattr(self, name, None)
            if widget is not None:
                widgets.append((key, widget))

        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            widgets.append(("lai_class", self.mMapLayerComboBox_LAI_Class))

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
            "version": 3,
            "mhm_version": (
                self.comboBox_mHMversion.currentText().strip()
                if hasattr(self, "comboBox_mHMversion") else ""
            ),
            "layers": layers,
            "text_inputs": text_inputs,
            "grid_resolutions": grid_resolutions,
            "grid_configuration": self.grid_configuration_snapshot(),
            "lai_input_type": (
                self.comboBox_laiInputType.currentText()
                if hasattr(self, "comboBox_laiInputType") else ""
            ),
            "folder_search": self.checkBox_enableFolderSearch.isChecked(),
            "categorical_types": {
                kind: self.categorical_input_mode(kind)
                for kind in ("lc", "soil", "geology")
            },
            "categorical_lookups": self._serialized_categorical_lookups(),
            "advanced_inputs": {
                kind: value.as_dict() if hasattr(value, "as_dict") else value
                for kind, value in self._advanced_inputs.items()
            },
            "land_cover_ready_source": self._land_cover_ready_source,
            "meteo_inputs": self.serialized_meteo_inputs(),
            "pour_point_outlet_id_field": self.selected_outlet_id_field(),
            "dem_domain": self.checkBox_DEMdomain.isChecked(),
            "domain_definition_type": (
                self.comboBox_domainDefinitionType.currentText().strip()
                if hasattr(self, "comboBox_domainDefinitionType") else ""
            ),
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
            path = saved.get("lookup_table")
            if path and self.project_folder:
                try:
                    relative = os.path.relpath(path, self.project_folder)
                    if relative != ".." and not relative.startswith(f"..{os.sep}"):
                        saved["lookup_table"] = relative.replace("\\", "/")
                except ValueError:
                    pass
            configs[kind] = saved
        return configs

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
        for kind in ("lc", "soil", "geology"):
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
        if hasattr(self, "comboBox_domainDefinitionType"):
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
        self._categorical_modes = {}
        self._advanced_inputs = {}
        self._land_cover_ready_source = ""
        self._domain_definition_mode = ""
        self._preferred_l1_resolution = None
        self._preferred_l11_resolution = None

    def restore_domain_inputs(self, state):
        """Restore the selected outlet ID field and DEM-domain flag."""
        self.populate_pour_point_outlet_fields(
            preferred=str(state.get("pour_point_outlet_id_field", "") or "")
        )
        self.checkBox_DEMdomain.blockSignals(True)
        self.checkBox_DEMdomain.setChecked(bool(state.get("dem_domain", False)))
        self.checkBox_DEMdomain.blockSignals(False)
        combo = getattr(self, "comboBox_domainDefinitionType", None)
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
            if kind not in {"lc", "soil", "geology"} or not isinstance(config, dict):
                continue
            restored = dict(config)
            path = restored.get("lookup_table")
            if path and not os.path.isabs(path):
                restored["lookup_table"] = os.path.abspath(
                    os.path.join(self.project_folder, path)
                )
            configs[kind] = restored
        self._categorical_lookup_configs = configs

        for kind, text in state.get("categorical_types", {}).items():
            if kind not in {"lc", "soil", "geology"}:
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
        combo_box = getattr(self, "comboBox_mHMversion", None)
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
            ("DEM", self.mMapLayerComboBox_dem),
            ("Pour points", self.mMapLayerComboBox_pour_points),
            ("Land cover", self.mMapLayerComboBox_land_cover),
            ("Soil", self.mMapLayerComboBox_soil),
            ("Geology", self.mMapLayerComboBox_geology),
        ]

        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            layer_inputs.append(("LAI class", self.mMapLayerComboBox_LAI_Class))

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

        if hasattr(self, "comboBox_laiInputType"):
            self.log_message(
                f"LAI input type: {self.comboBox_laiInputType.currentText() or '<not selected>'}"
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
            (getattr(self, "tab_geometry", None), getattr(self, "page_geometry", None)),
            (getattr(self, "tab_meteo", None), getattr(self, "page_meteo", None)),
            (getattr(self, "tab_hydro", None), getattr(self, "page_hydro", None)),
            (
                getattr(self, "tab_configuration", None),
                getattr(self, "page_configuration", None),
            ),
            (
                getattr(self, "tab_execution", None),
                getattr(self, "page_execution", None),
            ),
            (
                getattr(self, "tab_calibration", None),
                getattr(self, "page_calibration", None),
            ),
            (getattr(self, "tab_outputs", None), getattr(self, "page_outputs", None)),
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

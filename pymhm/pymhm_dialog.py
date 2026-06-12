# -*- coding: utf-8 -*-
"""Dialog and UI wiring for the pymhm QGIS plugin."""

import os
import json

# QGIS and PyQt imports
from qgis.PyQt.QtWidgets import QDialog, QFileDialog
from qgis.core import (
    QgsApplication,
    QgsMapLayer,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
)

# UI class from the compiled .ui file
from .ui_pymhm_dialog_base import Ui_pymhmDialog
from .qgis_compat import map_layer_filters

# Import utility mixin and processors
from .utils import DialogUtils
from .Morphology import MorphologyProcessor
from .Meteorology import MeteorologyProcessor
from .configuration_processor import ConfigurationProcessor
from .project_layout import (
    data_folder,
    data_raw_folder,
    ensure_project_structure,
    geometry_folder,
    output_folder,
    raw_meteo_folder,
    restart_folder,
    z_temp_folder,
)


class pymhmDialog(QDialog, Ui_pymhmDialog, DialogUtils):
    def __init__(self, parent=None):
        """Constructor."""

        super(pymhmDialog, self).__init__(parent)
        self.setupUi(self)

        # --- Filter map layer combo boxes to show only relevant layer types ---
        self.mMapLayerComboBox_dem.setFilters(map_layer_filters("RasterLayer"))
        self.mMapLayerComboBox_pour_points.setFilters(map_layer_filters("VectorLayer"))

        # Set filters for new layer combo boxes (both vector and raster allowed)
        self.mMapLayerComboBox_soil.setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.mMapLayerComboBox_land_cover.setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        self.mMapLayerComboBox_geology.setFilters(
            map_layer_filters("RasterLayer", "VectorLayer")
        )
        if hasattr(self, "mMapLayerComboBox_landCoverLookup"):
            self.mMapLayerComboBox_landCoverLookup.setFilters(
                map_layer_filters("VectorLayer")
            )
        for lookup_combo_name in (
            "mMapLayerComboBox_soilLookup",
            "mMapLayerComboBox_geologyLookup",
            "mMapLayerComboBox_laiLookup",
        ):
            lookup_combo = getattr(self, lookup_combo_name, None)
            if lookup_combo is not None:
                lookup_combo.setFilters(map_layer_filters("VectorLayer"))
        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            self.mMapLayerComboBox_LAI_Class.setFilters(
                map_layer_filters("RasterLayer", "VectorLayer")
            )
        self.configure_input_layer_combo_boxes()
        self.configure_lookup_field_combo_boxes()
        self.update_lai_input_controls()

        # --- Instance attributes for managing file paths ---
        self.project_folder = None
        self.geometry_folder = None  # Subfolder for geometry outputs
        self.input_state_filename = "pymhm_input_state.json"
        self._loading_input_state = False

        # --- Initialize processors ---
        self.morphology_processor = MorphologyProcessor(self)
        self.meteorology_processor = MeteorologyProcessor(self)
        self.configuration_processor = ConfigurationProcessor(self)
        self.simulation_processor = self.configuration_processor

        # --- Connect signals and slots ---
        self.configure_page_aliases()
        self.connect_signals()
        self.configuration_processor.refresh_status_indicators()

    def configure_page_aliases(self):
        """Keep renamed stacked pages compatible with older plugin code."""

        if hasattr(self, "page_configuration") and not hasattr(self, "page_validation"):
            self.page_validation = self.page_configuration
        if hasattr(self, "page_execution") and not hasattr(self, "page_datasets"):
            self.page_datasets = self.page_execution

    def configure_input_layer_combo_boxes(self):
        """Allow input layer boxes to start empty so layers are chosen deliberately."""

        layer_combo_boxes = [
            self.mMapLayerComboBox_dem,
            self.mMapLayerComboBox_pour_points,
            self.mMapLayerComboBox_soil,
            self.mMapLayerComboBox_land_cover,
            self.mMapLayerComboBox_geology,
        ]

        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            layer_combo_boxes.append(self.mMapLayerComboBox_LAI_Class)
        if hasattr(self, "mMapLayerComboBox_landCoverLookup"):
            layer_combo_boxes.append(self.mMapLayerComboBox_landCoverLookup)
        for combo_name in (
            "mMapLayerComboBox_soilLookup",
            "mMapLayerComboBox_geologyLookup",
            "mMapLayerComboBox_laiLookup",
        ):
            combo_box = getattr(self, combo_name, None)
            if combo_box is not None:
                layer_combo_boxes.append(combo_box)

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

    def connect_signals(self):
        """Connect all UI element signals to appropriate slots."""

        # Project management
        self.pushButton_BrowseProjectFolder.clicked.connect(self.select_project_folder)
        self.tabWidget.currentChanged.connect(self.on_tab_changed)

        # Morphology/Geometry processing - delegate to processor
        self.connect_processor_button(
            self.pushButton_convertDEMtoASC,
            "Convert DEM to ASC",
            self.morphology_processor.convert_dem_to_asc,
        )
        self.connect_processor_button(
            self.pushButton_fillDem, "Fill DEM", self.morphology_processor.fill_dem
        )
        self.connect_processor_button(
            self.pushButton_createNetwork,
            "Create Channel Network",
            self.morphology_processor.process_channel_network,
        )
        self.connect_processor_button(
            self.pushButton_snapPoints,
            "Snap Pour Points",
            self.morphology_processor.snap_points,
        )
        self.connect_processor_button(
            self.pushButton_delineate,
            "Delineate Watershed",
            self.morphology_processor.delineate_watershed,
        )
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
        self.connect_processor_button(
            self.pushButton_aspect, "Aspect", self.morphology_processor.process_aspect
        )
        self.connect_processor_button(
            self.pushButton_slope, "Slope", self.morphology_processor.process_slope
        )
        self.connect_processor_button(
            self.pushButton_flowAccumulation,
            "Flow Accumulation",
            self.morphology_processor.process_flow_accumulation,
        )
        self.connect_processor_button(
            self.pushButton_flowDirection,
            "Flow Direction",
            self.morphology_processor.process_flow_direction,
        )
        self.connect_processor_button(
            self.pushButton_gaugePosition,
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
                lambda layer=None: self.morphology_processor.update_gauged_outlet_count(
                    layer
                )
            )
        except Exception:
            pass
        self.morphology_processor.update_gauged_outlet_count()

        # Layer processing - delegate to processor
        self.connect_processor_button(
            self.pushButton_landUse,
            "Land Use",
            self.morphology_processor.process_land_use,
        )
        self.connect_processor_button(
            self.pushButton_soil, "Soil", self.morphology_processor.process_soil
        )
        self.connect_processor_button(
            self.pushButton_hydrogeology,
            "Hydrogeology",
            self.morphology_processor.process_geology,
        )

        # Mask all layers
        self.connect_processor_button(
            self.pushButton_maskAll,
            "Mask All Layers",
            self.morphology_processor.mask_all_layers,
        )

        # Write all layers (convert to ASCII)
        self.connect_processor_button(
            self.pushButton_writeAll,
            "Write All Layers",
            self.morphology_processor.write_all_layers,
        )

        # Reset geometry
        self.pushButton_resetGeometry.clicked.connect(
            self.morphology_processor.resetGeometry
        )

        # Execute all processing
        self.connect_processor_button(
            self.pushButton_executeAll,
            "Execute All Processing",
            self.morphology_processor.execute_all_processing,
        )

        # LAI file browser
        if hasattr(self, "pushButton_browse_lai"):
            self.pushButton_browse_lai.clicked.connect(self.browse_lai_file)

        # Lookup table browsers
        if hasattr(self, "pushButton_browse_land_cover_lookup"):
            self.pushButton_browse_land_cover_lookup.clicked.connect(
                self.browse_land_cover_lookup
            )
        if hasattr(self, "pushButton_browse_soil_lookup"):
            self.pushButton_browse_soil_lookup.clicked.connect(self.browse_soil_lookup)
        if hasattr(self, "pushButton_browse_geology_lookup"):
            self.pushButton_browse_geology_lookup.clicked.connect(
                self.browse_geology_lookup
            )
        self.connect_lookup_field_signals()

        # Meteorology folder browser
        self.pushButton_browse_meteo_folder.clicked.connect(self.browse_meteo_folder)
        if hasattr(self, "pushButton_meteo_runSave"):
            self.connect_processor_button(
                self.pushButton_meteo_runSave,
                "Prepare Meteorology Forcing",
                self.meteorology_processor.process_meteo_forcing,
            )

        # Configuration/execution processing - delegate to processor
        if hasattr(self, "pushButton_configureMHM"):
            self.connect_processor_button(
                self.pushButton_configureMHM,
                "Configure mHM",
                self.configuration_processor.configure_mhm,
            )
        if hasattr(self, "pushButton_configureParameters"):
            self.connect_processor_button(
                self.pushButton_configureParameters,
                "Configure Parameters",
                self.configuration_processor.configure_parameters,
            )
        if hasattr(self, "pushButton_configureOutputs"):
            self.connect_processor_button(
                self.pushButton_configureOutputs,
                "Configure Outputs",
                self.configuration_processor.configure_outputs,
            )
        self.connect_processor_button(
            self.pushButton_createNML,
            "Create Namelists",
            self.configuration_processor.create_nml_files,
        )
        self.connect_processor_button(
            self.pushButton_RUN, "Run mHM", self.configuration_processor.run_mhm
        )
        if hasattr(self, "comboBox_mHMversion"):
            self.comboBox_mHMversion.currentIndexChanged.connect(
                lambda index=None: self.configuration_processor.handle_version_changed()
            )
        if hasattr(self, "pushButton_browseConfiguration"):
            self.pushButton_browseConfiguration.clicked.connect(
                self.configuration_processor.browse_configuration_file
            )
        self.connect_input_state_signals()

        # Initialize CRS widget with project CRS
        project_crs = QgsProject.instance().crs()
        if project_crs.isValid():
            self.mProjectionSelectionWidget_crs.setCrs(project_crs)

    def connect_input_state_signals(self):
        """Save input selections when editable inputs change."""

        for _, combo_box in self.input_layer_widgets():
            try:
                combo_box.layerChanged.connect(
                    lambda layer=None: self.save_input_state()
                )
            except Exception:
                pass

        for _, line_edit in self.input_text_widgets():
            try:
                line_edit.editingFinished.connect(self.save_input_state)
            except Exception:
                try:
                    line_edit.textChanged.connect(
                        lambda text=None: self.save_input_state()
                    )
                except Exception:
                    pass

        for _, combo_box in self.input_lookup_field_widgets():
            try:
                combo_box.currentTextChanged.connect(
                    lambda text=None: self.save_input_state()
                )
            except Exception:
                try:
                    combo_box.currentIndexChanged.connect(
                        lambda index=None: self.save_input_state()
                    )
                except Exception:
                    pass

        if hasattr(self, "comboBox_laiInputType"):
            try:
                self.comboBox_laiInputType.currentIndexChanged.connect(
                    lambda index=None: self.save_input_state()
                )
            except Exception:
                pass

        try:
            self.mProjectionSelectionWidget_crs.crsChanged.connect(
                lambda crs=None: self.save_input_state()
            )
        except Exception:
            pass

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

    def lookup_field_specs(self):
        """Return lookup table layer widgets and their companion field widgets."""
        return [
            (
                "land_cover_lookup_field",
                "mMapLayerComboBox_landCoverLookup",
                "comboBox_landCoverLookupField",
            ),
            (
                "soil_lookup_field",
                "mMapLayerComboBox_soilLookup",
                "comboBox_soilLookupField",
            ),
            (
                "geology_lookup_field",
                "mMapLayerComboBox_geologyLookup",
                "comboBox_geologyLookupField",
            ),
            (
                "lai_lookup_field",
                "mMapLayerComboBox_laiLookup",
                "comboBox_laiLookupField",
            ),
        ]

    def input_lookup_field_widgets(self):
        """Return lookup field combo boxes and their state keys."""
        widgets = []
        for key, _, field_combo_name in self.lookup_field_specs():
            field_combo = getattr(self, field_combo_name, None)
            if field_combo is not None:
                widgets.append((key, field_combo))
        return widgets

    def configure_lookup_field_combo_boxes(self):
        """Populate lookup field dropdowns from their selected table layers."""
        for _, layer_combo_name, field_combo_name in self.lookup_field_specs():
            layer_combo = getattr(self, layer_combo_name, None)
            field_combo = getattr(self, field_combo_name, None)
            if layer_combo is None or field_combo is None:
                continue
            self.populate_lookup_field_combo(
                field_combo,
                layer_combo.currentLayer(),
                field_combo.currentText(),
            )

    def connect_lookup_field_signals(self):
        """Refresh lookup field dropdowns when lookup table layers change."""
        for _, layer_combo_name, field_combo_name in self.lookup_field_specs():
            layer_combo = getattr(self, layer_combo_name, None)
            field_combo = getattr(self, field_combo_name, None)
            if layer_combo is None or field_combo is None:
                continue
            try:
                layer_combo.layerChanged.connect(
                    lambda layer=None, combo=field_combo: self.populate_lookup_field_combo(
                        combo, layer
                    )
                )
            except Exception:
                pass

        if hasattr(self, "comboBox_laiInputType"):
            try:
                self.comboBox_laiInputType.currentIndexChanged.connect(
                    lambda index=None: self.update_lai_input_controls()
                )
            except Exception:
                pass

    def populate_lookup_field_combo(self, field_combo, lookup_layer=None, preferred_field=""):
        """Load a lookup field combo box with fields from a QGIS table layer."""
        current_field = preferred_field or field_combo.currentText()
        try:
            field_combo.blockSignals(True)
        except Exception:
            pass

        field_combo.clear()
        field_names = []
        if lookup_layer is not None and isinstance(lookup_layer, QgsVectorLayer):
            try:
                field_names = lookup_layer.fields().names()
            except Exception:
                field_names = []

        field_combo.addItems(field_names)
        if current_field and current_field in field_names:
            field_combo.setCurrentText(current_field)
        elif field_names:
            field_combo.setCurrentIndex(0)

        field_combo.setEnabled(bool(field_names))
        try:
            field_combo.blockSignals(False)
        except Exception:
            pass
        self.update_lai_input_controls()

    def restore_lookup_fields(self, lookup_fields):
        """Restore saved lookup table field selections after layers are restored."""
        for key, layer_combo_name, field_combo_name in self.lookup_field_specs():
            layer_combo = getattr(self, layer_combo_name, None)
            field_combo = getattr(self, field_combo_name, None)
            if layer_combo is None or field_combo is None:
                continue
            self.populate_lookup_field_combo(
                field_combo,
                layer_combo.currentLayer(),
                lookup_fields.get(key, ""),
            )

    def restore_lai_input_type(self, lai_input_type):
        """Restore the selected LAI input type by text when possible."""
        combo_box = getattr(self, "comboBox_laiInputType", None)
        if combo_box is None or not lai_input_type:
            return
        index = combo_box.findText(lai_input_type)
        if index >= 0:
            combo_box.setCurrentIndex(index)

    def lai_uses_lookup_table(self):
        """Return True when the LAI workflow uses class and lookup-table inputs."""
        combo_box = getattr(self, "comboBox_laiInputType", None)
        if combo_box is None:
            return False
        return combo_box.currentText().strip().lower() == "lai classes and lookup table"

    def update_lai_input_controls(self):
        """Enable LAI lookup controls only for the class-and-lookup workflow."""
        enabled = self.lai_uses_lookup_table()
        lookup_combo = getattr(self, "mMapLayerComboBox_laiLookup", None)
        field_combo = getattr(self, "comboBox_laiLookupField", None)
        if lookup_combo is not None:
            lookup_combo.setEnabled(enabled)
        if field_combo is not None:
            field_combo.setEnabled(
                enabled and field_combo.count() > 0
            )

    def input_layer_widgets(self):
        """Return persistent layer input widgets and their state keys."""
        widgets = [
            ("dem", self.mMapLayerComboBox_dem),
            ("pour_points", self.mMapLayerComboBox_pour_points),
            ("soil", self.mMapLayerComboBox_soil),
            ("land_cover", self.mMapLayerComboBox_land_cover),
            ("geology", self.mMapLayerComboBox_geology),
        ]

        if hasattr(self, "mMapLayerComboBox_LAI_Class"):
            widgets.append(("lai_class", self.mMapLayerComboBox_LAI_Class))
        if hasattr(self, "mMapLayerComboBox_landCoverLookup"):
            widgets.append(
                ("land_cover_lookup", self.mMapLayerComboBox_landCoverLookup)
            )
        for key, combo_name in (
            ("soil_lookup", "mMapLayerComboBox_soilLookup"),
            ("geology_lookup", "mMapLayerComboBox_geologyLookup"),
            ("lai_lookup", "mMapLayerComboBox_laiLookup"),
        ):
            combo_box = getattr(self, combo_name, None)
            if combo_box is not None:
                widgets.append((key, combo_box))

        return widgets

    def input_text_widgets(self):
        """Return persistent text/path input widgets and their state keys."""
        widgets = []

        for key, widget_name in (
            ("lai_file", "lineEdit_lai_file"),
            ("soil_lookup_file", "lineEdit_soil_lookup"),
            ("geology_lookup_file", "lineEdit_geology_lookup"),
            ("meteo_folder", "lineEdit_meteo_folder"),
        ):
            widget = getattr(self, widget_name, None)
            if widget is not None:
                widgets.append((key, widget))

        if hasattr(self, "lineEdit_land_cover_lookup"):
            widgets.append(("land_cover_lookup_file", self.lineEdit_land_cover_lookup))
        if hasattr(self, "lineEdit_loadConfiguration"):
            widgets.append(
                ("configuration_settings_file", self.lineEdit_loadConfiguration)
            )

        for key, widget_name in (
            ("l0", "lineEdit_L0"),
            ("l1", "lineEdit_L1"),
            ("l2", "lineEdit_L2"),
        ):
            if hasattr(self, widget_name):
                widgets.append((key, getattr(self, widget_name)))

        return widgets

    def input_state_path(self):
        """Return the project-local input state file path."""
        if not self.project_folder:
            return None
        return os.path.join(self.project_folder, self.input_state_filename)

    def save_input_state(self):
        """Save selected inputs to a JSON file in the project folder."""
        if not self.project_folder or self._loading_input_state:
            return

        state_path = self.input_state_path()
        if not state_path:
            return

        layers = {}
        for key, combo_box in self.input_layer_widgets():
            layer = combo_box.currentLayer()
            if not layer:
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
            }

        text_inputs = {key: widget.text() for key, widget in self.input_text_widgets()}
        lookup_fields = {
            key: combo_box.currentText()
            for key, combo_box in self.input_lookup_field_widgets()
        }

        crs = self.get_crs()
        state = {
            "version": 1,
            "layers": layers,
            "text_inputs": text_inputs,
            "lookup_fields": lookup_fields,
            "lai_input_type": (
                self.comboBox_laiInputType.currentText()
                if hasattr(self, "comboBox_laiInputType") else ""
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
            with open(state_path, "w", encoding="utf-8") as state_file:
                json.dump(state, state_file, indent=2, sort_keys=True)
        except Exception as e:
            self.log_message(f"WARNING: Could not save input state: {e}")

    def load_input_state(self):
        """Load saved input selections from the project folder."""
        state_path = self.input_state_path()
        if not state_path or not os.path.exists(state_path):
            self.log_message("No saved input state found for this project.")
            return

        try:
            with open(state_path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
        except Exception as e:
            self.log_message(f"WARNING: Could not read input state: {e}")
            return

        self._loading_input_state = True
        try:
            self.restore_text_inputs(state.get("text_inputs", {}))
            self.restore_lai_input_type(state.get("lai_input_type", ""))
            self.restore_input_crs(state.get("crs_authid", ""))
            self.restore_input_layers(state.get("layers", {}))
            self.restore_lookup_fields(state.get("lookup_fields", {}))
            self.update_lai_input_controls()
        finally:
            self._loading_input_state = False

        self.log_message(f"Input state loaded: {state_path}")

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
        if hasattr(self, "mMapLayerComboBox_landCoverLookup"):
            layer_inputs.append(
                ("Land cover lookup", self.mMapLayerComboBox_landCoverLookup)
            )
        for label, combo_name in (
            ("Soil lookup", "mMapLayerComboBox_soilLookup"),
            ("Geology lookup", "mMapLayerComboBox_geologyLookup"),
            ("LAI lookup", "mMapLayerComboBox_laiLookup"),
        ):
            combo_box = getattr(self, combo_name, None)
            if combo_box is not None:
                layer_inputs.append((label, combo_box))

        for label, combo_box in layer_inputs:
            layer = combo_box.currentLayer()
            if layer:
                self.log_message(f"{label}: {layer.name()} | {layer.source()}")
            else:
                self.log_message(f"{label}: <not selected>")

        if hasattr(self, "comboBox_laiInputType"):
            self.log_message(
                f"LAI input type: {self.comboBox_laiInputType.currentText() or '<not selected>'}"
            )
        if hasattr(self, "lineEdit_lai_file"):
            self.log_message(
                f"LAI file: {self.lineEdit_lai_file.text() or '<not selected>'}"
            )
        for label, combo_name in (
            ("Land cover lookup field", "comboBox_landCoverLookupField"),
            ("Soil lookup field", "comboBox_soilLookupField"),
            ("Geology lookup field", "comboBox_geologyLookupField"),
            ("LAI lookup field", "comboBox_laiLookupField"),
        ):
            combo_box = getattr(self, combo_name, None)
            if combo_box is not None and combo_box.isEnabled():
                self.log_message(
                    f"{label}: {combo_box.currentText() or '<not selected>'}"
                )
        if hasattr(self, "lineEdit_land_cover_lookup"):
            self.log_message(
                f"Land cover lookup: {self.lineEdit_land_cover_lookup.text() or '<not selected>'}"
            )
        if hasattr(self, "lineEdit_soil_lookup"):
            self.log_message(
                f"Soil lookup: {self.lineEdit_soil_lookup.text() or '<not selected>'}"
            )
        if hasattr(self, "lineEdit_geology_lookup"):
            self.log_message(
                f"Geology lookup: {self.lineEdit_geology_lookup.text() or '<not selected>'}"
            )
        if hasattr(self, "lineEdit_meteo_folder"):
            self.log_message(
                f"Meteorology folder: {self.lineEdit_meteo_folder.text() or '<not selected>'}"
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

            self.load_input_state()
            if not self.lineEdit_meteo_folder.text().strip():
                default_meteo_folder = raw_meteo_folder(self.project_folder)
                self.lineEdit_meteo_folder.setText(default_meteo_folder)
                self.log_message(
                    f"Meteorology data folder set to: {default_meteo_folder}"
                )
            self.morphology_processor.update_gauged_outlet_count()

            # Load project state in morphology processor
            self.morphology_processor.load_project_state()
            self.meteorology_processor.load_project_state()
            self.configuration_processor.load_project_state()
            self.configuration_processor.refresh_status_indicators()

    def on_tab_changed(self, index):
        """Switches the stacked widget page when the tab is changed."""
        page = self.page_for_tab_index(index)
        if page is not None:
            self.stackedWidget.setCurrentWidget(page)
        else:
            self.stackedWidget.setCurrentIndex(index)
        self.configuration_processor.refresh_status_indicators()
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

    def browse_land_cover_lookup(self):
        """Browse for Land Cover lookup table or database"""
        if not hasattr(self, "lineEdit_land_cover_lookup"):
            self.log_message(
                "Land cover lookup is now selected from the layer dropdown."
            )
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Land Cover Lookup Table",
            "",
            "All Files (*);;CSV (*.csv);;Excel (*.xlsx *.xls);;Database (*.db *.sqlite)",
        )
        if file_path:
            self.lineEdit_land_cover_lookup.setText(file_path)
            self.log_message(f"Land cover lookup table selected: {file_path}")
            self.save_input_state()

    def browse_soil_lookup(self):
        """Browse for Soil lookup table or database"""
        if not hasattr(self, "lineEdit_soil_lookup"):
            self.log_message("Soil lookup is now selected from the layer dropdown.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Soil Lookup Table",
            "",
            "All Files (*);;CSV (*.csv);;Excel (*.xlsx *.xls);;Database (*.db *.sqlite)",
        )
        if file_path:
            self.lineEdit_soil_lookup.setText(file_path)
            self.log_message(f"Soil lookup table selected: {file_path}")
            self.save_input_state()

    def browse_geology_lookup(self):
        """Browse for Geology lookup table or database"""
        if not hasattr(self, "lineEdit_geology_lookup"):
            self.log_message("Geology lookup is now selected from the layer dropdown.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Geology Lookup Table",
            "",
            "All Files (*);;CSV (*.csv);;Excel (*.xlsx *.xls);;Database (*.db *.sqlite)",
        )
        if file_path:
            self.lineEdit_geology_lookup.setText(file_path)
            self.log_message(f"Geology lookup table selected: {file_path}")
            self.save_input_state()

    def browse_meteo_folder(self):
        """Browse for Meteorology data folder"""
        folder = QFileDialog.getExistingDirectory(
            self, "Select Meteorology Data Folder"
        )
        if folder:
            self.lineEdit_meteo_folder.setText(folder)
            self.log_message(f"Meteorology data folder selected: {folder}")
            self.save_input_state()

    def get_lai_time_range(self):
        """Get the selected LAI time/date for extraction"""
        if hasattr(self, "dateEdit") and self.dateEdit:
            return self.dateEdit.dateTime()
        return None

    def get_crs(self):
        """Get the selected CRS"""
        return self.mProjectionSelectionWidget_crs.crs()

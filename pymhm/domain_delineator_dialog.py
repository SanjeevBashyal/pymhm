# -*- coding: utf-8 -*-
"""Interactive per-outlet watershed delineation dialog."""
from __future__ import annotations

import os
import re
from pathlib import Path

from qgis.PyQt.QtCore import QEvent, Qt
from qgis.PyQt.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QMessageBox,
    QVBoxLayout,
)
from qgis.core import (
    QgsCoordinateTransform,
    QgsGeometry,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
)
from qgis.gui import QgsMapCanvas, QgsMapToolEmitPoint

from .input_selection import scan_project_inputs
from .Morphology.hydrology.discharge_writer import (
    write_streamflow_observation,
)
from .Morphology.hydrology.observation_paths import (
    streamflow_observation_folder,
)
from .Morphology.hydrology.outlets import (
    StationIdError,
    station_id_int,
    station_id_text,
)
from .Morphology.watershed.domain_state import (
    load_state,
    resolve_input_path,
    resolve_output_path,
    save_state,
)
from .project_layout import geometry_folder
from .pyui.ui_domain_delineator_dialog import Ui_DomainDelineatorDialog


class DomainDelineatorDialog(QDialog, Ui_DomainDelineatorDialog):
    """Configure gauges and delineate one optional domain per outlet."""

    def __init__(
        self,
        main_dialog,
        processor,
        pour_points_layer,
        outlet_id_field,
        outlet_ids,
    ):
        super().__init__(main_dialog)
        self.setupUi(self)
        self.main_dialog = main_dialog
        self.processor = processor
        self.project_folder = str(main_dialog.project_folder)
        self.pour_points_layer = pour_points_layer
        self.outlet_id_field = str(outlet_id_field)
        self.outlet_ids = [str(value) for value in outlet_ids]
        self.current_outlet_id = ""
        self._preview_result = None
        self._preview_channel_path = ""
        self._map_tool = None
        self._picking = False

        self.canvas = QgsMapCanvas(self.widget_mapCanvasHost)
        layout = QVBoxLayout(self.widget_mapCanvasHost)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)

        self._load_state()
        self._features = self._outlet_features()
        self._prepare_map_layers()
        self._connect_signals()
        self.finished.connect(self._cleanup)
        self.listWidget_outlets.addItems(self.outlet_ids)
        if self.outlet_ids:
            self.listWidget_outlets.setCurrentRow(0)

    def _connect_signals(self):
        self.listWidget_outlets.currentTextChanged.connect(
            self._load_outlet)
        self.checkBox_isGaugedOutlet.toggled.connect(
            self._set_discharge_enabled)
        self.checkBox_isDomainOutlet.toggled.connect(
            self.widget_domainControls.setEnabled)
        self.pushButton_browseDischargeFile.clicked.connect(
            self._browse_discharge)
        self.pushButton_generateChannelNetwork.clicked.connect(
            self._generate_channel_network)
        self.pushButton_pickLocation.clicked.connect(self._start_picking)
        self.pushButton_save.clicked.connect(self._save_outlet)
        self.pushButton_close.clicked.connect(self.reject)

    def _load_state(self):
        source = self._layer_source(self.pour_points_layer)
        state = load_state(self.project_folder)
        if (
            state.get("pour_points_source")
            and state["pour_points_source"] != source
        ) or (
            state.get("outlet_id_field")
            and state["outlet_id_field"] != self.outlet_id_field
        ):
            state["outlets"] = {}

        records = state.setdefault("outlets", {})
        state["outlets"] = {
            outlet_id: dict(records.get(outlet_id, {}))
            for outlet_id in self.outlet_ids
        }
        state["pour_points_source"] = source
        state["outlet_id_field"] = self.outlet_id_field
        state["dem_domain"] = bool(self.main_dialog.checkBox_DEMdomain.isChecked())
        self.state = state

    @staticmethod
    def _layer_source(layer):
        source = str(layer.source() or "")
        local_path = source.split("|", 1)[0]
        if os.path.exists(local_path):
            return str(Path(local_path).resolve())
        return source

    def _outlet_features(self):
        features = {}
        for feature in self.pour_points_layer.getFeatures():
            outlet_id = station_id_text(feature.attribute(self.outlet_id_field))
            if outlet_id:
                features[outlet_id] = feature
        return features

    def _prepare_map_layers(self):
        if not self.processor._ensure_filled_dem(self.processor.fill_dem):
            raise RuntimeError("The filled DEM could not be prepared.")
        if not self.processor._ensure_flow_accumulation(
            self.processor.process_flow_accumulation,
            self.processor.fill_dem,
        ):
            raise RuntimeError("Flow accumulation could not be prepared.")
        if not self.processor._ensure_channel_network(
            self.processor.process_channel_network,
            self.processor.process_flow_accumulation,
            self.processor.fill_dem,
        ):
            raise RuntimeError("The channel network could not be prepared.")

        self._filled_dem_layer = QgsRasterLayer(
            self.processor.filled_dem_path, "Filled DEM")
        self._flow_accumulation_layer = QgsRasterLayer(
            self.processor.flow_accumulation_path, "Flow accumulation")
        self._channel_layer = QgsVectorLayer(
            self.processor.channel_network_vector_path,
            "Channel network",
            "ogr",
        )
        if not self._filled_dem_layer.isValid():
            raise RuntimeError("The prepared filled DEM is invalid.")
        if not self._flow_accumulation_layer.isValid():
            raise RuntimeError("The prepared flow accumulation is invalid.")
        if not self._channel_layer.isValid():
            raise RuntimeError("The prepared channel network is invalid.")
        self._flow_context = self.processor._build_flwdir_from_filled_dem()
        if not self._flow_context:
            raise RuntimeError("The filled DEM flow grid could not be prepared.")
        self.canvas.setDestinationCrs(self._filled_dem_layer.crs())
        self._watershed_layer = None
        self._refresh_canvas()

    def _refresh_canvas(self):
        layers = []
        for layer in (
            self._watershed_layer,
            self.pour_points_layer,
            self._channel_layer,
            self._flow_accumulation_layer,
        ):
            if layer is not None and layer.isValid():
                layers.append(layer)
        self.canvas.setLayers(layers)
        self.canvas.refresh()

    def _load_outlet(self, outlet_id):
        self.current_outlet_id = str(outlet_id or "")
        if not self.current_outlet_id:
            return

        self._preview_result = None
        record = self.state["outlets"].get(self.current_outlet_id, {})
        self._channel_layer = QgsVectorLayer(
            self.processor.channel_network_vector_path,
            "Channel network",
            "ogr",
        )
        self.label_outletIDValue.setText(self.current_outlet_id)
        self.checkBox_isGaugedOutlet.blockSignals(True)
        self.checkBox_isDomainOutlet.blockSignals(True)
        self.checkBox_isGaugedOutlet.setChecked(
            bool(record.get("is_gauged", record.get("gauged", False))))
        self.checkBox_isDomainOutlet.setChecked(
            bool(record.get("is_domain", record.get("domain", False))))
        self.checkBox_isGaugedOutlet.blockSignals(False)
        self.checkBox_isDomainOutlet.blockSignals(False)
        self._set_discharge_enabled(
            self.checkBox_isGaugedOutlet.isChecked())
        self.widget_domainControls.setEnabled(
            self.checkBox_isDomainOutlet.isChecked())
        try:
            threshold = max(1, int(record.get("threshold_cells", 1) or 1))
        except (TypeError, ValueError):
            threshold = 1
        self.spinBox_channelThreshold.setValue(threshold)

        discharge = self._resolved_input(record.get("discharge_file"))
        self._populate_discharge_files(discharge)
        area = record.get("catchment_area_m2")
        self._show_area(area)
        self._show_saved_watershed(record)
        self._zoom_to_outlet()

    def _set_discharge_enabled(self, enabled):
        self.comboBox_dischargeFile.setEnabled(bool(enabled))
        self.pushButton_browseDischargeFile.setEnabled(bool(enabled))

    @staticmethod
    def _normal_path(path):
        return os.path.normcase(os.path.abspath(str(path)))

    def _resolved_input(self, value):
        if not value:
            return ""
        return str(resolve_input_path(self.project_folder, value))

    def _used_discharge_paths(self, exclude_outlet=None):
        paths = set(self.main_dialog.selected_input_file_paths())
        for outlet_id, record in self.state["outlets"].items():
            if outlet_id == exclude_outlet or not record.get("discharge_file"):
                continue
            paths.add(self._normal_path(self._resolved_input(
                record["discharge_file"])))
        return paths

    def _populate_discharge_files(self, current_path=""):
        current = self._normal_path(current_path) if current_path else ""
        excluded = self._used_discharge_paths(self.current_outlet_id)
        self.comboBox_dischargeFile.blockSignals(True)
        self.comboBox_dischargeFile.clear()
        self.comboBox_dischargeFile.addItem("", "")

        available = {}
        for item in scan_project_inputs(self.project_folder, "lookup"):
            path = self._normal_path(item.data["path"])
            if path in excluded:
                continue
            available[path] = item.label
        if current and current not in excluded and os.path.isfile(current):
            available.setdefault(current, self._path_label(current))

        selected_index = 0
        for path, label in sorted(
            available.items(), key=lambda item: item[1].casefold()
        ):
            self.comboBox_dischargeFile.addItem(label, path)
            if path == current:
                selected_index = self.comboBox_dischargeFile.count() - 1
        self.comboBox_dischargeFile.setCurrentIndex(selected_index)
        self.comboBox_dischargeFile.blockSignals(False)

    def _path_label(self, path):
        try:
            return Path(path).resolve().relative_to(
                Path(self.project_folder).resolve()).as_posix()
        except ValueError:
            return str(Path(path).resolve())

    def _browse_discharge(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select discharge data",
            self.project_folder,
            "Discharge data (*.csv *.txt)",
        )
        if not path:
            return
        if Path(path).suffix.lower() not in {".csv", ".txt"}:
            QMessageBox.warning(
                self, "Invalid file", "Select a CSV or TXT discharge file.")
            return

        normalized = self._normal_path(path)
        if normalized in self._used_discharge_paths(self.current_outlet_id):
            QMessageBox.warning(
                self,
                "File already selected",
                "This file is already selected elsewhere in the plugin.",
            )
            return
        self._select_discharge_path(normalized)

    def _select_discharge_path(self, path):
        for index in range(self.comboBox_dischargeFile.count()):
            if self.comboBox_dischargeFile.itemData(index) == path:
                self.comboBox_dischargeFile.setCurrentIndex(index)
                return
        self.comboBox_dischargeFile.addItem(self._path_label(path), path)
        self.comboBox_dischargeFile.setCurrentIndex(
            self.comboBox_dischargeFile.count() - 1)

    def _generate_channel_network(self):
        output_folder = self._domain_output_folder()
        output_path = os.path.join(
            output_folder, "2_channel_network_preview.shp")
        try:
            self._channel_layer = QgsVectorLayer(
                self.processor.channel_network_vector_path,
                "Channel network",
                "ogr",
            )
            self._refresh_canvas()
            QApplication.setOverrideCursor(Qt.WaitCursor)
            result = self.processor.create_channel_network(
                self.spinBox_channelThreshold.value(),
                output_path,
                force=True,
                load=False,
            )
            if not result:
                raise RuntimeError("Channel network generation failed.")
            self._preview_channel_path = result
            self._channel_layer = QgsVectorLayer(
                result, "Channel network preview", "ogr")
            if not self._channel_layer.isValid():
                raise RuntimeError("The generated channel network is invalid.")
            self._refresh_canvas()
        except Exception as error:
            self._show_error("Channel Network", error)
        finally:
            QApplication.restoreOverrideCursor()

    def _start_picking(self):
        if not self.current_outlet_id:
            return
        self._stop_picking()
        self._map_tool = QgsMapToolEmitPoint(self.canvas)
        self._map_tool.canvasClicked.connect(self._point_picked)
        self.canvas.setMapTool(self._map_tool)
        self.canvas.setCursor(Qt.CrossCursor)
        self.pushButton_pickLocation.setText("Click on map…")
        self._picking = True
        QApplication.instance().installEventFilter(self)

    def eventFilter(self, watched, event):
        if self._picking and event.type() == QEvent.MouseButtonPress:
            position = self.canvas.mapFromGlobal(event.globalPos())
            if not self.canvas.rect().contains(position):
                self._stop_picking()
                QMessageBox.warning(
                    self,
                    "Pick location",
                    "Click inside the map canvas to select an outlet location.",
                )
                return True
        return super().eventFilter(watched, event)

    def _stop_picking(self):
        if not self._picking and self._map_tool is None:
            return
        self._picking = False
        application = QApplication.instance()
        if application is not None:
            application.removeEventFilter(self)
        if self._map_tool is not None:
            try:
                self._map_tool.canvasClicked.disconnect(self._point_picked)
            except (TypeError, RuntimeError):
                pass
            self.canvas.unsetMapTool(self._map_tool)
        self._map_tool = None
        self.canvas.unsetCursor()
        self.pushButton_pickLocation.setText("Pick location")

    def _point_picked(self, point, button):
        if button != Qt.LeftButton:
            return
        self._stop_picking()
        raster_path, vector_path = self._outlet_paths(
            self.current_outlet_id, preview=True)
        try:
            self._watershed_layer = None
            self._refresh_canvas()
            QApplication.setOverrideCursor(Qt.WaitCursor)
            result = self.processor.delineate_single_outlet(
                point.x(),
                point.y(),
                raster_path,
                vector_path,
                basin_id=self.outlet_ids.index(self.current_outlet_id) + 1,
                context=self._flow_context,
            )
            if not result:
                raise RuntimeError("Watershed delineation failed.")
            center_x, center_y = result["cell_center"]
            result["picked"] = {
                "x": float(center_x),
                "y": float(center_y),
                "crs": self._filled_dem_layer.crs().authid(),
            }
            self._preview_result = result
            self._show_area(result["catchment_area_m2"])
            self._watershed_layer = QgsVectorLayer(
                result["vector_path"],
                f"Watershed {self.current_outlet_id}",
                "ogr",
            )
            self._refresh_canvas()
        except Exception as error:
            self._show_error("Watershed Delineation", error)
        finally:
            QApplication.restoreOverrideCursor()

    def _save_outlet(self):
        outlet_id = self.current_outlet_id
        if not outlet_id:
            return
        record = dict(self.state["outlets"].get(outlet_id, {}))
        was_gauged = bool(
            record.get("is_gauged", record.get("gauged", False)))
        is_gauged = self.checkBox_isGaugedOutlet.isChecked()
        is_domain = self.checkBox_isDomainOutlet.isChecked()
        discharge_path = str(
            self.comboBox_dischargeFile.currentData() or "")
        cursor_set = False

        try:
            if is_gauged:
                station_id_int(outlet_id)
                if not discharge_path:
                    raise ValueError(
                        "Select a discharge CSV or TXT file for this gauge.")
                if self._normal_path(discharge_path) in (
                    self._used_discharge_paths(outlet_id)
                ):
                    raise ValueError(
                        "The selected discharge file is already used elsewhere.")

            QApplication.setOverrideCursor(Qt.WaitCursor)
            cursor_set = True
            result = None
            if is_domain:
                picked = (
                    self._preview_result.get("picked")
                    if self._preview_result
                    else record.get("picked")
                )
                if not isinstance(picked, dict):
                    raise ValueError(
                        "Pick a map location before saving this domain.")
                x, y = self._picked_in_dem_crs(picked)
                raster_path, vector_path = self._outlet_paths(outlet_id)
                self._watershed_layer = None
                self._refresh_canvas()
                result = self.processor.delineate_single_outlet(
                    x,
                    y,
                    raster_path,
                    vector_path,
                    basin_id=self.outlet_ids.index(outlet_id) + 1,
                    context=self._flow_context,
                )
                if not result:
                    raise RuntimeError("Watershed delineation failed.")

            record["is_gauged"] = is_gauged
            record["is_domain"] = is_domain
            record["threshold_cells"] = self.spinBox_channelThreshold.value()
            record["discharge_file"] = discharge_path if is_gauged else ""
            if result:
                center_x, center_y = result["cell_center"]
                record.update({
                    "picked": {
                        "x": float(center_x),
                        "y": float(center_y),
                        "crs": self._filled_dem_layer.crs().authid(),
                    },
                    "catchment_area_m2": result["catchment_area_m2"],
                    "mask_path": result["raster_path"],
                    "vector_path": result["vector_path"],
                })

            self.state["outlets"][outlet_id] = record
            self.state["dem_domain"] = bool(
                self.main_dialog.checkBox_DEMdomain.isChecked())
            if self.state["dem_domain"]:
                self._prepare_dem_domain()
            self._merge_active_domains()

            if is_gauged:
                layer = QgsVectorLayer(
                    discharge_path, Path(discharge_path).name, "ogr")
                output = write_streamflow_observation(
                    layer,
                    outlet_id,
                    Path(streamflow_observation_folder(self.project_folder)),
                )
                self.processor.mark_output_prepared(
                    str(output), name=output.name, loaded=False)

            save_state(self.project_folder, self.state)
            if was_gauged and not is_gauged:
                self._remove_observation_file(outlet_id)
            self.main_dialog.save_input_state()
            self.processor.update_gauged_outlet_count()
            self._preview_result = None
            self._show_saved_watershed(record)
            QApplication.restoreOverrideCursor()
            cursor_set = False
            QMessageBox.information(
                self, "Domain Delineator", f"Saved outlet {outlet_id}.")
        except (StationIdError, ValueError, RuntimeError) as error:
            self._show_error("Save Outlet", error)
        except Exception as error:
            self._show_error("Save Outlet", error)
        finally:
            if cursor_set:
                QApplication.restoreOverrideCursor()

    def _remove_observation_file(self, outlet_id):
        path = (
            Path(streamflow_observation_folder(self.project_folder))
            / f"{outlet_id}.txt"
        )
        try:
            path.unlink(missing_ok=True)
            key = self.processor.output_state_key(str(path))
            self.processor.processing_state.get("outputs", {}).pop(key, None)
            self.processor.save_processing_state()
        except OSError as error:
            self.main_dialog.log_message(
                f"WARNING: Could not remove old discharge output: {error}")

    def _picked_in_dem_crs(self, picked):
        x = float(picked["x"])
        y = float(picked["y"])
        authid = str(picked.get("crs", "") or "")
        target = self._filled_dem_layer.crs()
        if not authid or authid == target.authid():
            return x, y

        from qgis.core import QgsCoordinateReferenceSystem, QgsPointXY

        source = QgsCoordinateReferenceSystem(authid)
        if not source.isValid():
            raise ValueError("The saved outlet coordinate CRS is invalid.")
        transform = QgsCoordinateTransform(
            source, target, QgsProject.instance())
        point = transform.transform(QgsPointXY(x, y))
        return point.x(), point.y()

    def _prepare_dem_domain(self):
        output_raster, output_vector = self._dem_paths()
        reference = self.processor._read_raster_array(
            self.processor.filled_dem_path, as_float=True)
        if not reference:
            raise RuntimeError("Could not read the filled DEM.")
        deps = self.processor._get_python_morphology_deps()
        dem, invalid_mask, _ = self.processor._normalise_dem_array(
            reference["array"], reference["nodata"])
        if dem is None:
            raise RuntimeError("Could not determine the valid DEM domain.")
        domain = deps["np"].where(~invalid_mask, 1, 0).astype(
            deps["np"].int32)
        if not self.processor._write_raster_array(
            output_raster,
            domain,
            reference,
            nodata=0,
            gdal_type=deps["gdal"].GDT_Int32,
        ):
            raise RuntimeError("Could not write the DEM-domain mask.")
        if not self._polygonize_mask(output_raster, output_vector):
            raise RuntimeError("Could not write the DEM-domain polygon.")

    def _polygonize_mask(self, raster_path, vector_path):
        raw_path = os.path.splitext(vector_path)[0] + "_raw.shp"
        self.processor._remove_vector_dataset(raw_path)
        try:
            result = self.processor.run_processing_algorithm(
                "gdal:polygonize",
                {
                    "INPUT": raster_path,
                    "BAND": 1,
                    "FIELD": "DN",
                    "EIGHT_CONNECTEDNESS": False,
                    "EXTRA": "",
                    "OUTPUT": raw_path,
                },
            )
            return bool(
                result
                and os.path.exists(raw_path)
                and self.processor._copy_nonzero_polygons(
                    raw_path, vector_path)
            )
        finally:
            self.processor._remove_vector_dataset(raw_path)

    def _merge_active_domains(self):
        vector_paths = []
        for outlet_id, record in self.state["outlets"].items():
            if not record.get("is_domain", record.get("domain", False)):
                continue
            value = record.get("vector_path")
            path = str(resolve_output_path(
                self.project_folder, value)) if value else ""
            if not path or not os.path.exists(path):
                raise ValueError(
                    f"Delineate and save domain outlet {outlet_id} first.")
            vector_paths.append(path)

        if self.state.get("dem_domain"):
            dem_vector = self._dem_paths()[1]
            if not os.path.exists(dem_vector):
                raise RuntimeError("The DEM-domain polygon is missing.")
            vector_paths.append(dem_vector)

        merged_path = os.path.join(
            geometry_folder(self.project_folder),
            "Watersheds",
            "4_watershed_merged_vector.shp",
        )
        if not vector_paths:
            self.processor._remove_vector_dataset(merged_path)
            self.processor.merged_watershed_path = None
            self.processor.watershed_vector_path = None
            return

        for path in vector_paths:
            layer = QgsVectorLayer(path, Path(path).stem, "ogr")
            if not layer.isValid() or layer.featureCount() < 1:
                raise ValueError(
                    f"Domain polygon is invalid or empty: {Path(path).name}")

        pending_path = os.path.splitext(merged_path)[0] + "_pending.shp"
        self.processor._remove_vector_dataset(pending_path)
        try:
            result = self.main_dialog.run_processing_algorithm(
                "native:mergevectorlayers",
                {"LAYERS": vector_paths, "CRS": None, "OUTPUT": pending_path},
            )
            pending_layer = QgsVectorLayer(
                pending_path, "Pending merged domains", "ogr")
            valid = (
                bool(result)
                and pending_layer.isValid()
                and pending_layer.featureCount() > 0
            )
            pending_layer = None
            if not valid:
                raise RuntimeError(
                    "Could not merge the active domain polygons.")
            self._publish_vector_dataset(pending_path, merged_path)
        finally:
            self.processor._remove_vector_dataset(pending_path)

        self.processor.merged_watershed_path = merged_path
        self.processor.watershed_vector_path = merged_path
        self.processor.mark_output_prepared(
            merged_path, name="4_watershed_merged", loaded=False)

    def _publish_vector_dataset(self, source_path, target_path):
        source_base = os.path.splitext(source_path)[0]
        target_base = os.path.splitext(target_path)[0]
        self.processor._remove_vector_dataset(target_path)
        for extension in (
            ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj", ".fix"
        ):
            source = source_base + extension
            if os.path.exists(source):
                os.replace(source, target_base + extension)

    def _domain_output_folder(self):
        path = os.path.join(
            geometry_folder(self.project_folder),
            "Watersheds",
            "DomainDelineation",
        )
        os.makedirs(path, exist_ok=True)
        return path

    def _safe_outlet_name(self, outlet_id):
        value = re.sub(r"[^A-Za-z0-9_-]+", "_", str(outlet_id)).strip("_")
        index = self.outlet_ids.index(str(outlet_id)) + 1
        return f"{index}_{value or 'outlet'}"

    def _outlet_paths(self, outlet_id, preview=False):
        prefix = "_preview_" if preview else "4_watershed_"
        base = os.path.join(
            self._domain_output_folder(),
            prefix + self._safe_outlet_name(outlet_id),
        )
        return base + ".tif", base + ".shp"

    def _dem_paths(self):
        base = os.path.join(
            self._domain_output_folder(), "4_watershed_DEM")
        return base + ".tif", base + ".shp"

    def _show_saved_watershed(self, record):
        self._watershed_layer = None
        value = record.get("vector_path")
        if value:
            try:
                path = str(resolve_output_path(self.project_folder, value))
                layer = QgsVectorLayer(
                    path, f"Watershed {self.current_outlet_id}", "ogr")
                if layer.isValid():
                    self._watershed_layer = layer
            except ValueError:
                pass
        self._refresh_canvas()

    def _show_area(self, value):
        try:
            area = float(value)
        except (TypeError, ValueError):
            self.label_catchmentAreaValue.setText("-")
            return
        self.label_catchmentAreaValue.setText(
            f"{area / 1_000_000.0:.3f} km²")

    def _zoom_to_outlet(self):
        feature = self._features.get(self.current_outlet_id)
        if feature is None or feature.geometry().isEmpty():
            return
        geometry = QgsGeometry(feature.geometry())
        source = self.pour_points_layer.crs()
        target = self._filled_dem_layer.crs()
        if source.isValid() and target.isValid() and source != target:
            transform = QgsCoordinateTransform(
                source, target, QgsProject.instance())
            geometry.transform(transform)

        extent = geometry.boundingBox()
        padding = max(
            self._filled_dem_layer.extent().width(),
            self._filled_dem_layer.extent().height(),
        ) / 40.0
        extent.grow(padding if padding > 0 else 1.0)
        self.canvas.setExtent(extent)
        self.canvas.refresh()

    def _show_error(self, title, error):
        self.main_dialog.log_message(f"ERROR: {title}: {error}")
        QMessageBox.critical(self, title, str(error))

    def _cleanup(self, *_args):
        self._stop_picking()
        self.canvas.setLayers([])

    def closeEvent(self, event):
        self._cleanup()
        super().closeEvent(event)

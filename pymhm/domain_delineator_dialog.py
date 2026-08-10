# -*- coding: utf-8 -*-
"""Interactive per-outlet watershed delineation dialog."""
from __future__ import annotations

import copy
import os
from pathlib import Path

from qgis.PyQt.QtCore import QEvent, Qt
from qgis.PyQt.QtGui import QColor
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
from qgis.gui import QgsMapCanvas, QgsMapToolEmitPoint, QgsVertexMarker

from .input_selection import scan_project_inputs
from .Morphology.hydrology.discharge_dialog import OutletAssignment
from .Morphology.hydrology.outlets import (
    StationIdError,
    station_id_text,
)
from .Morphology.watershed.domain_state import (
    DOMAIN_MODE_DELINEATOR,
    active_domain_records,
    assign_domain_ids,
    resolve_input_path,
    resolve_output_path,
    save_state,
)
from .Morphology.watershed.domain_workflow import DomainWorkflow
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
        self.workflow = DomainWorkflow(
            main_dialog,
            processor,
            pour_points_layer,
            outlet_id_field,
            self.outlet_ids,
        )
        self.current_outlet_id = ""
        self._preview_result = None
        self._preview_channel_path = ""
        self._map_tool = None
        self._picking = False

        self.canvas = QgsMapCanvas(self.widget_mapCanvasHost)
        layout = QVBoxLayout(self.widget_mapCanvasHost)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self._outlet_marker = QgsVertexMarker(self.canvas)
        self._outlet_marker.setIconType(QgsVertexMarker.ICON_CIRCLE)
        self._outlet_marker.setIconSize(22)
        self._outlet_marker.setPenWidth(4)
        self._outlet_marker.setColor(QColor(220, 30, 30))
        self._outlet_marker.setFillColor(QColor(255, 240, 40, 180))
        self._outlet_marker.hide()

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
        self.state = self.workflow.load_synced_state(
            DOMAIN_MODE_DELINEATOR,
            bool(self.main_dialog.checkBox_DEMdomain.isChecked()),
        )

    @staticmethod
    def _layer_source(layer):
        return DomainWorkflow.layer_source(layer)

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
            self._outlet_marker.hide()
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
        proposed_state = copy.deepcopy(self.state)
        record = dict(proposed_state["outlets"].get(outlet_id, {}))
        was_gauged = bool(
            record.get("is_gauged", record.get("gauged", False)))
        is_gauged = self.checkBox_isGaugedOutlet.isChecked()
        is_domain = self.checkBox_isDomainOutlet.isChecked()
        discharge_path = str(
            self.comboBox_dischargeFile.currentData() or "")
        cursor_set = False

        try:
            assignment = OutletAssignment(
                outlet_id=outlet_id,
                is_gauge=is_gauged,
                is_domain=is_domain,
                discharge_layer=(
                    QgsVectorLayer(
                        discharge_path,
                        Path(discharge_path).name,
                        "ogr",
                    )
                    if is_gauged and discharge_path
                    else None
                ),
            )
            if is_gauged:
                if not discharge_path:
                    raise ValueError(
                        "Select a discharge CSV or TXT file for this gauge.")
                if self._normal_path(discharge_path) in (
                    self._used_discharge_paths(outlet_id)
                ):
                    raise ValueError(
                        "The selected discharge file is already used elsewhere.")
            prepared = self.workflow.validate_gauge_assignments([assignment])

            self.workflow.apply_assignment_records(
                proposed_state,
                [assignment],
                prepared,
            )
            record = proposed_state["outlets"][outlet_id]
            record["threshold_cells"] = self.spinBox_channelThreshold.value()
            proposed_state["dem_domain"] = bool(
                self.main_dialog.checkBox_DEMdomain.isChecked()
            )
            assign_domain_ids(proposed_state)
            self.workflow.require_active_domain(proposed_state)
            self.workflow.validate_unique_state_gauge_ids(proposed_state)

            if is_domain:
                picked = (
                    self._preview_result.get("picked")
                    if self._preview_result
                    else record.get("picked")
                )
                if not isinstance(picked, dict):
                    raise ValueError(
                        "Pick a map location before saving this domain.")
                record["picked"] = dict(picked)

            delineations = []
            for domain in active_domain_records(proposed_state):
                if domain.get("is_dem_domain"):
                    continue
                domain_outlet_id = domain["outlet_id"]
                domain_record = proposed_state["outlets"][domain_outlet_id]
                picked = domain_record.get("picked")
                if not isinstance(picked, dict):
                    raise ValueError(
                        "Pick and save a location for domain outlet "
                        f"{domain_outlet_id}."
                    )
                x, y = self._picked_in_dem_crs(picked)
                raster_path, vector_path = self._outlet_paths(
                    domain_outlet_id,
                    state=proposed_state,
                )
                delineations.append(
                    (
                        domain_outlet_id,
                        domain_record,
                        x,
                        y,
                        raster_path,
                        vector_path,
                    )
                )

            QApplication.setOverrideCursor(Qt.WaitCursor)
            cursor_set = True
            self._watershed_layer = None
            self._refresh_canvas()
            for (
                domain_outlet_id,
                domain_record,
                x,
                y,
                raster_path,
                vector_path,
            ) in delineations:
                result = self.processor.delineate_single_outlet(
                    x,
                    y,
                    raster_path,
                    vector_path,
                    basin_id=int(domain_record["domain_id"]),
                    context=self._flow_context,
                )
                if not result:
                    raise RuntimeError(
                        "Watershed delineation failed for outlet "
                        f"{domain_outlet_id}."
                    )
                center_x, center_y = result["cell_center"]
                domain_record.update(
                    {
                        "picked": {
                            "x": float(center_x),
                            "y": float(center_y),
                            "crs": self._filled_dem_layer.crs().authid(),
                        },
                        "catchment_area_m2": result["catchment_area_m2"],
                        "mask_path": result["raster_path"],
                        "vector_path": result["vector_path"],
                    }
                )

            if proposed_state["dem_domain"]:
                self.workflow.prepare_dem_domain(proposed_state)
            self.workflow.merge_active_domains(proposed_state)
            self.workflow.update_gauge_domain_ids(
                proposed_state,
                self.pour_points_layer,
            )
            self.workflow.write_gauges(prepared)
            save_state(self.project_folder, proposed_state)
            self.state = proposed_state
            if was_gauged and not is_gauged:
                self._remove_observation_file(outlet_id)
            self.main_dialog.save_input_state()
            self.processor.update_gauged_outlet_count()
            self._preview_result = None
            self._show_saved_watershed(record if is_domain else {})
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
        self.workflow.remove_observation_file(outlet_id)

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
        return self.workflow.prepare_dem_domain(self.state)

    def _merge_active_domains(self):
        return self.workflow.merge_active_domains(self.state)

    def _domain_output_folder(self):
        return self.workflow.domain_output_folder()

    def _outlet_paths(self, outlet_id, preview=False, state=None):
        return self.workflow.outlet_paths(
            outlet_id,
            state or self.state,
            preview=preview,
        )

    def _dem_paths(self):
        return self.workflow.dem_paths()

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
            self._outlet_marker.hide()
            return
        geometry = QgsGeometry(feature.geometry())
        source = self.pour_points_layer.crs()
        target = self._filled_dem_layer.crs()
        if source.isValid() and target.isValid() and source != target:
            transform = QgsCoordinateTransform(
                source, target, QgsProject.instance())
            geometry.transform(transform)

        marker_point = geometry.asPoint()
        if marker_point.isEmpty():
            marker_point = geometry.centroid().asPoint()
        self._outlet_marker.setCenter(marker_point)
        self._outlet_marker.show()

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
        if self._outlet_marker is not None:
            self._outlet_marker.hide()
            try:
                self.canvas.scene().removeItem(self._outlet_marker)
            except RuntimeError:
                pass
            self._outlet_marker = None
        self.canvas.setLayers([])

    def closeEvent(self, event):
        self._cleanup()
        super().closeEvent(event)

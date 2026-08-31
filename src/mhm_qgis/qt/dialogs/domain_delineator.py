# -*- coding: utf-8 -*-
"""Interactive per-outlet watershed delineation dialog."""
from __future__ import annotations

import copy
import os
from functools import partial
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
from ...qt.dialogs.discharge_assignment import OutletAssignment
from ...Morphology.hydrology.outlets import (
    StationIdError,
    station_id_text,
)
from ...core.handlers.state.domain_state import (
    DOMAIN_MODE_DELINEATOR,
    active_domain_records,
    assign_domain_ids,
    resolve_input_path,
    resolve_output_path,
    save_state,
)
from ...Morphology.watershed.domain_workflow import DomainWorkflow
from ...Morphology.file_tasks import (
    delineate_domains_file,
    delineate_outlet_file,
)
from ...project_layout import geometry_folder
from ...project_layout import domain_data_folder, domain_dem_path
from ...qt.controllers import domain_delineator as domain_delineator_controller
from ...qt.bindings.domain_delineator import bind as bind_domain_delineator
from ...qt.ui.pyui.ui_domain_delineator_dialog import Ui_DomainDelineatorDialog
from ...viewport_raster_range import ViewportRasterRangeController


def _run_outlet_task(task, options):
    return delineate_outlet_file(task=task, **options)


def _run_domains_task(task, options):
    return delineate_domains_file(task=task, **options)


class DomainDelineatorDialog(QDialog, Ui_DomainDelineatorDialog):
    """Configure gauges and delineate one optional domain per outlet."""

    def __init__(
        self,
        main_dialog,
        processor,
        pour_points_layer,
        outlet_id_field,
        outlet_ids,
        prepared_context=None,
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
        self.canvas.setEnabled(True)
        self.widget_mapCanvasHost.setEnabled(True)
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
        self._draft_state = copy.deepcopy(self.state)
        self._features = self._outlet_features()
        self._prepare_map_layers(prepared_context)
        self._connect_signals()
        self.finished.connect(self._cleanup)
        self.listWidget_outlets.addItems(self.outlet_ids)
        if self.outlet_ids:
            self.listWidget_outlets.setCurrentRow(0)

    def _connect_signals(self):
        bind_domain_delineator(self)

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

    def _prepare_map_layers(self, prepared_context=None):
        if prepared_context is None:
            raise RuntimeError("Domain Delineator prerequisites were not prepared.")
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
        self._flow_context = None
        self.canvas.setDestinationCrs(self._filled_dem_layer.crs())
        self._viewport_range = ViewportRasterRangeController(
            self.canvas,
            self._flow_accumulation_layer,
            self.processor.flow_accumulation_path,
            self.main_dialog.task_coordinator,
            self,
        )
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

    def _load_outlet(self, *args, **kwargs):
        return domain_delineator_controller._load_outlet(self, *args, **kwargs)

    def _set_discharge_enabled(self, *args, **kwargs):
        return domain_delineator_controller._set_discharge_enabled(self, *args, **kwargs)

    @staticmethod
    def _normal_path(path):
        return os.path.normcase(os.path.abspath(str(path)))

    def _resolved_input(self, value):
        if not value:
            return ""
        return str(resolve_input_path(self.project_folder, value))

    def _used_discharge_paths(self, exclude_outlet=None):
        paths = set(self.main_dialog.selected_input_file_paths())
        for outlet_id, record in self._draft_state["outlets"].items():
            if outlet_id == exclude_outlet or not record.get("discharge_file"):
                continue
            paths.add(self._normal_path(self._resolved_input(
                record["discharge_file"])))
        return paths

    def _populate_discharge_files(self, *args, **kwargs):
        return domain_delineator_controller._populate_discharge_files(self, *args, **kwargs)

    def _path_label(self, *args, **kwargs):
        return domain_delineator_controller._path_label(self, *args, **kwargs)

    def _browse_discharge(self, *args, **kwargs):
        return domain_delineator_controller._browse_discharge(self, *args, **kwargs)

    def _select_discharge_path(self, *args, **kwargs):
        return domain_delineator_controller._select_discharge_path(self, *args, **kwargs)

    def _generate_channel_network(self):
        output_folder = self._domain_output_folder()
        output_path = os.path.join(
            output_folder, "2_channel_network_preview.shp")
        self._channel_layer = QgsVectorLayer(
            self.processor.channel_network_vector_path,
            "Channel network",
            "ogr",
        )
        self._refresh_canvas()
        threshold = self.spinBox_channelThreshold.value()

        self.main_dialog.morphology_tasks.start_hydrology(
            key="domain-channel-preview",
            controls=(
                self.pushButton_generateChannelNetwork,
                self.pushButton_close,
            ),
            load="",
            channel_path=output_path,
            threshold=threshold,
            direction=False,
            write_flow=False,
            done=lambda result: self._channel_network_ready(
                result["channel_network"]
            ),
            failed=lambda message: self._show_error(
                "Channel Network", str(message).split("\n", 1)[0]
            ),
        )

    def _channel_network_ready(self, result):
        self._preview_channel_path = str(result)
        self._channel_layer = QgsVectorLayer(
            str(result), "Channel network preview", "ogr"
        )
        if not self._channel_layer.isValid():
            self._show_error("Channel Network", "The generated channel network is invalid.")
            return
        self._refresh_canvas()

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
        self._watershed_layer = None
        self._refresh_canvas()
        outlet_id = self.current_outlet_id
        x, y = point.x(), point.y()

        options = {
            "filled_dem": self.processor.filled_dem_path,
            "x": x,
            "y": y,
            "raster_path": raster_path,
            "vector_path": vector_path,
            "basin_id": self.outlet_ids.index(outlet_id) + 1,
        }
        self.main_dialog.task_coordinator.submit(
            "domain-watershed-preview",
            f"Preview watershed {outlet_id}",
            partial(_run_outlet_task, options=options),
            controls=(
                self.pushButton_pickLocation,
                self.pushButton_nextPourPoint,
                self.pushButton_save,
                self.listWidget_outlets,
                self.pushButton_close,
            ),
            on_success=lambda result: self._watershed_preview_ready(
                outlet_id, result
            ),
            on_error=lambda message: self._show_error(
                "Watershed Delineation", str(message).split("\n", 1)[0]
            ),
            resource="morphology-files",
            task_aware=True,
        )

    def _watershed_preview_ready(self, outlet_id, result):
        if outlet_id != self.current_outlet_id:
            return
        try:
            center_x, center_y = result["cell_center"]
            result["picked"] = {
                "x": float(center_x),
                "y": float(center_y),
                "crs": self._filled_dem_layer.crs().authid(),
            }
            self._preview_result = result
            self._show_area(result["catchment_area_m2"])
            self._show_picked_coordinates(result["picked"])
            self._watershed_layer = QgsVectorLayer(
                result["vector_path"],
                f"Watershed {outlet_id}",
                "ogr",
            )
            self._refresh_canvas()
        except Exception as error:
            self._show_error("Watershed Delineation", error)

    def _current_assignment(self, *args, **kwargs):
        return domain_delineator_controller._current_assignment(self, *args, **kwargs)

    def _stage_current_outlet(self):
        """Validate and retain one outlet without publishing final outputs."""
        if not self.current_outlet_id:
            return
        assignment = self._current_assignment()
        if assignment.is_gauge:
            path = self._normal_path(
                self.comboBox_dischargeFile.currentData() or ""
            )
            if path in self._used_discharge_paths(self.current_outlet_id):
                raise ValueError("The selected discharge file is already used elsewhere.")
        prepared = self.workflow.validate_gauge_assignments([assignment])
        self.workflow.apply_assignment_records(
            self._draft_state, [assignment], prepared
        )
        record = self._draft_state["outlets"][self.current_outlet_id]
        record["threshold_cells"] = self.spinBox_channelThreshold.value()
        picked = (
            self._preview_result.get("picked")
            if self._preview_result
            else record.get("picked")
        )
        if (assignment.is_domain or assignment.is_gauge) and not isinstance(
            picked, dict
        ):
            picked = self._default_gauge_point(self.current_outlet_id)
        if isinstance(picked, dict):
            record["picked"] = dict(picked)
            if assignment.is_gauge:
                record["gauge_point"] = dict(picked)
                record["gauge_point"]["source"] = "picked"
        elif assignment.is_gauge:
            record["gauge_point"] = self._default_gauge_point(
                self.current_outlet_id
            )
        else:
            record.pop("gauge_point", None)
        self._draft_state["dem_domain"] = bool(
            self.main_dialog.checkBox_DEMdomain.isChecked()
        )
        assign_domain_ids(self._draft_state)
        self._preview_result = None

    def _default_gauge_point(self, outlet_id):
        feature = self._features.get(str(outlet_id))
        if feature is None or feature.geometry().isEmpty():
            raise ValueError(f"Gauge outlet {outlet_id} has no default point location.")
        geometry = QgsGeometry(feature.geometry())
        source = self.pour_points_layer.crs()
        target = self._filled_dem_layer.crs()
        if source.isValid() and target.isValid() and source != target:
            transform = QgsCoordinateTransform(
                source, target, QgsProject.instance()
            )
            transform.setBallparkTransformsAreAppropriate(True)
            geometry.transform(transform)
        point = geometry.asPoint()
        if point.isEmpty():
            point = geometry.centroid().asPoint()
        return {
            "x": float(point.x()),
            "y": float(point.y()),
            "crs": target.authid() if target.isValid() else "",
            "source": "default",
        }

    def _next_outlet(self, *args, **kwargs):
        return domain_delineator_controller._next_outlet(self, *args, **kwargs)

    def _all_assignments(self):
        assignments = []
        for outlet_id in self.outlet_ids:
            record = self._draft_state["outlets"].get(outlet_id, {})
            discharge_path = self._resolved_input(record.get("discharge_file"))
            assignments.append(
                OutletAssignment(
                    outlet_id=outlet_id,
                    is_gauge=bool(record.get("is_gauged", False)),
                    is_domain=bool(record.get("is_domain", False)),
                    discharge_layer=(
                        QgsVectorLayer(
                            discharge_path,
                            Path(discharge_path).name,
                            "ogr",
                        )
                        if record.get("is_gauged") and discharge_path
                        else None
                    ),
                )
            )
        return assignments

    def _save_outlet(self):
        try:
            self._stage_current_outlet()
            proposed_state = copy.deepcopy(self._draft_state)
            assignments = self._all_assignments()
            prepared = self.workflow.validate_gauge_assignments(assignments)
            previously_gauged = {
                outlet_id
                for outlet_id, record in self.state.get("outlets", {}).items()
                if record.get("is_gauged", False)
            }
            self.workflow.apply_assignment_records(
                proposed_state, assignments, prepared
            )
            assign_domain_ids(proposed_state)
            self.workflow.require_active_domain(proposed_state)
            self.workflow.validate_unique_state_gauge_ids(proposed_state)
            delineations = []
            for domain in active_domain_records(proposed_state):
                if domain.get("is_dem_domain"):
                    continue
                outlet_id = domain["outlet_id"]
                record = proposed_state["outlets"][outlet_id]
                picked = record.get("picked")
                if not isinstance(picked, dict):
                    raise ValueError(
                        f"Pick a map location for domain outlet {outlet_id}."
                    )
                raster_path, vector_path = self._outlet_paths(
                    outlet_id, state=proposed_state
                )
                x, y = self._picked_in_dem_crs(picked)
                delineations.append(
                    (
                        outlet_id,
                        x,
                        y,
                        raster_path,
                        vector_path,
                        int(record["domain_id"]),
                        domain_dem_path(self.project_folder, outlet_id),
                    )
                )
            dem_crs_authid = self._filled_dem_layer.crs().authid()
            dem_domain = None
            if proposed_state.get("dem_domain"):
                dem_raster, dem_vector = self.workflow.dem_paths()
                dem_domain = (
                    int(proposed_state["dem_domain_id"]),
                    dem_raster,
                    dem_vector,
                    domain_dem_path(self.project_folder, "dem_extent"),
                )
            merged_path = os.path.join(
                geometry_folder(self.project_folder),
                "Watersheds",
                "4_watershed_merged_vector.shp",
            )
            self._watershed_layer = None
            self._refresh_canvas()
            options = {
                "filled_dem": self.processor.filled_dem_path,
                "delineations": tuple(delineations),
                "dem_domain": dem_domain,
                "merged_path": merged_path,
            }
            started = self.main_dialog.task_coordinator.submit(
                "domain-final-save",
                "Save all domain and gauge inputs",
                partial(_run_domains_task, options=options),
                task_aware=True,
                controls=(
                    self.pushButton_nextPourPoint,
                    self.pushButton_save,
                    self.pushButton_pickLocation,
                    self.pushButton_generateChannelNetwork,
                    self.listWidget_outlets,
                    self.pushButton_close,
                ),
                on_success=lambda result: self._domain_files_ready(
                        result,
                        proposed_state,
                        prepared,
                        previously_gauged,
                        dem_crs_authid,
                ),
                on_error=lambda message: self._show_error(
                    "Save Domain Inputs", str(message).split("\n", 1)[0]
                ),
                resource="morphology-files",
            )
            if not started:
                raise RuntimeError("Domain saving is already running.")
        except (StationIdError, ValueError, RuntimeError) as error:
            self._show_error("Save Domain Inputs", error)
        except Exception as error:
            self._show_error("Save Domain Inputs", error)

    def _commit_domain_state(
        self,
        result,
        proposed_state,
        prepared,
        previously_gauged,
        dem_crs_authid,
    ):
        for outlet_id, delineation in result["outlets"].items():
            record = proposed_state["outlets"][outlet_id]
            center_x, center_y = delineation["cell_center"]
            record.update(
                {
                    "picked": {
                        "x": float(center_x),
                        "y": float(center_y),
                        "crs": dem_crs_authid,
                    },
                    "catchment_area_m2": delineation["catchment_area_m2"],
                    "mask_path": delineation["raster_path"],
                    "vector_path": delineation["vector_path"],
                    "domain_directory": domain_data_folder(
                        self.project_folder, outlet_id
                    ),
                    "dem_path": delineation.get("dem_path", ""),
                }
            )
        if result.get("dem_domain_path"):
            proposed_state["dem_domain_directory"] = domain_data_folder(
                self.project_folder, "dem_extent"
            )
            proposed_state["dem_domain_path"] = result["dem_domain_path"]
        self.processor.merged_watershed_path = result["merged_path"]
        self.processor.watershed_vector_path = result["merged_path"]
        self.processor.mark_output_prepared(
            result["merged_path"], name="4_watershed_merged", loaded=False
        )
        self.workflow.update_gauge_domain_ids(
            proposed_state, self.pour_points_layer
        )
        self.workflow.write_gauges(prepared)
        save_state(self.project_folder, proposed_state)
        self.workflow.remove_deselected_gauges(
            previously_gauged, proposed_state
        )
        from ...core.handlers.state.nml_settings import sync_domain_settings

        sync_domain_settings(self.project_folder)
        self._domain_state_saved(proposed_state)

    def _domain_files_ready(
        self,
        result,
        proposed_state,
        prepared,
        previously_gauged,
        dem_crs_authid,
    ):
        try:
            self._commit_domain_state(
                result,
                proposed_state,
                prepared,
                previously_gauged,
                dem_crs_authid,
            )
        except Exception as error:
            self._show_error("Save Domain Inputs", error)

    def _domain_state_saved(self, *args, **kwargs):
        return domain_delineator_controller._domain_state_saved(self, *args, **kwargs)

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

    def _show_area(self, *args, **kwargs):
        return domain_delineator_controller._show_area(self, *args, **kwargs)

    def _show_picked_coordinates(self, *args, **kwargs):
        return domain_delineator_controller._show_picked_coordinates(self, *args, **kwargs)

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

    def _show_error(self, *args, **kwargs):
        return domain_delineator_controller._show_error(self, *args, **kwargs)

    def _cleanup(self, *_args):
        self._stop_picking()
        controller = getattr(self, "_viewport_range", None)
        if controller is not None:
            controller.close()
            self._viewport_range = None
        if self._outlet_marker is not None:
            self._outlet_marker.hide()
            try:
                self.canvas.scene().removeItem(self._outlet_marker)
            except RuntimeError:
                pass
            self._outlet_marker = None
        self.canvas.setLayers([])

    def closeEvent(self, event):
        if self._domain_jobs_busy():
            event.ignore()
            QMessageBox.information(
                self,
                "Processing",
                "Wait for the active Domain Delineator task to finish.",
            )
            return
        self._cleanup()
        super().closeEvent(event)

    def reject(self):
        if self._domain_jobs_busy():
            QMessageBox.information(
                self,
                "Processing",
                "Wait for the active Domain Delineator task to finish.",
            )
            return
        super().reject()

    def _domain_jobs_busy(self):
        coordinator = self.main_dialog.task_coordinator
        return any(
            coordinator.is_busy(key)
            for key in (
                "domain-channel-preview",
                "domain-channel-preview-dem",
                "domain-watershed-preview",
                "domain-final-save",
            )
        )

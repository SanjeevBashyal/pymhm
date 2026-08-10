# -*- coding: utf-8 -*-
"""Shared domain and gauge assignment workflow."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from qgis.core import (
    QgsCoordinateTransform,
    QgsGeometry,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
)

from ...project_layout import geometry_folder
from ..hydrology.discharge_dialog import OutletAssignment
from ..hydrology.discharge_writer import (
    local_source_path,
    records_from_layer,
    streamflow_filename,
    write_streamflow_observation,
)
from ..hydrology.observation_paths import streamflow_observation_folder
from ..hydrology.outlets import (
    StationIdError,
    find_outlet_id_field,
    station_id_int,
    station_id_text,
)
from .domain_state import (
    DOMAIN_MODE_DEM_EXTENT,
    DOMAIN_MODE_SNAPPED,
    active_domain_records,
    assign_domain_ids,
    gauged_outlet_ids,
    load_state,
    resolve_output_path,
    save_state,
)


@dataclass(frozen=True)
class _PreparedGauge:
    assignment: OutletAssignment
    gauge_id: int
    input_path: Path
    output_path: Path


class DomainWorkflow:
    """Apply domain/gauge choices using an existing morphology processor."""

    def __init__(
        self,
        main_dialog,
        processor,
        pour_points_layer,
        outlet_id_field: str,
        outlet_ids: Iterable[str],
    ) -> None:
        self.main_dialog = main_dialog
        self.processor = processor
        self.project_folder = str(main_dialog.project_folder)
        self.pour_points_layer = pour_points_layer
        self.outlet_id_field = str(outlet_id_field)
        self.outlet_ids = [str(value) for value in outlet_ids]

    def load_synced_state(
        self,
        definition_mode: str = "",
        dem_domain: bool | None = None,
    ) -> dict:
        """Load state and align it with the current pour-point selection."""
        source = self.layer_source(self.pour_points_layer)
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
        state["outlet_order"] = list(self.outlet_ids)
        state["pour_points_source"] = source
        state["outlet_id_field"] = self.outlet_id_field
        if definition_mode:
            state["definition_mode"] = definition_mode
        if dem_domain is not None:
            state["dem_domain"] = bool(dem_domain)
        return assign_domain_ids(state)

    def apply_dem_extent(
        self,
        assignments: Iterable[OutletAssignment],
    ) -> dict:
        """Create one valid-DEM domain and commit optional gauges."""
        assignments = self._validated_assignment_list(assignments)
        prepared = self.validate_gauge_assignments(assignments)
        previously_gauged = set(
            gauged_outlet_ids(load_state(self.project_folder))
        )
        state = self.load_synced_state(DOMAIN_MODE_DEM_EXTENT, True)
        self.apply_assignment_records(state, assignments, prepared)
        for record in state["outlets"].values():
            record["is_domain"] = False
        assign_domain_ids(state)

        if not self.processor._ensure_filled_dem(self.processor.fill_dem):
            raise RuntimeError("The filled DEM could not be prepared.")
        self.prepare_dem_domain(state)
        self.merge_active_domains(state)
        self.update_gauge_domain_ids(state, self.pour_points_layer)
        self.write_gauges(prepared)
        save_state(self.project_folder, state)
        self.remove_deselected_gauges(previously_gauged, state)
        return state

    def apply_snapped_domains(
        self,
        assignments: Iterable[OutletAssignment],
        *,
        include_dem_domain: bool = False,
    ) -> dict:
        """Delineate only checked domains at their snapped pour points."""
        assignments = self._validated_assignment_list(assignments)
        prepared = self.validate_gauge_assignments(assignments)
        previously_gauged = set(
            gauged_outlet_ids(load_state(self.project_folder))
        )
        state = self.load_synced_state(
            DOMAIN_MODE_SNAPPED,
            include_dem_domain,
        )
        self.apply_assignment_records(state, assignments, prepared)
        assign_domain_ids(state)
        self.require_active_domain(state)

        self.regenerate_snapped_points()

        snapped_layer = QgsVectorLayer(
            self.processor.snapped_points_path,
            "Snapped pour points",
            "ogr",
        )
        if not snapped_layer.isValid():
            raise RuntimeError("The prepared snapped pour-point layer is invalid.")
        features = self._features_by_outlet(snapped_layer, require_exact=True)
        context = self.processor._build_flwdir_from_filled_dem()
        if not context:
            raise RuntimeError("The filled DEM flow grid could not be prepared.")
        transform = self.processor._snapped_to_filled_dem_transform(snapped_layer)
        status_field = self.processor._snap_status_field(
            snapped_layer.fields().names()
        )
        for assignment in assignments:
            if not (assignment.is_domain or assignment.is_gauge):
                continue
            feature = features[assignment.outlet_id]
            if status_field and str(feature.attribute(status_field)) == "failed":
                raise ValueError(
                    f"Outlet {assignment.outlet_id} could not be snapped "
                    "to the channel network."
                )

        for assignment in assignments:
            record = state["outlets"][assignment.outlet_id]
            if not assignment.is_domain:
                continue
            feature = features[assignment.outlet_id]
            point = feature.geometry().asPoint()
            if transform is not None:
                point = transform.transform(point)
            raster_path, vector_path = self.outlet_paths(
                assignment.outlet_id,
                state,
            )
            result = self.processor.delineate_single_outlet(
                point.x(),
                point.y(),
                raster_path,
                vector_path,
                basin_id=int(record["domain_id"]),
                context=context,
            )
            if not result:
                raise RuntimeError(
                    f"Watershed delineation failed for outlet {assignment.outlet_id}."
                )
            center_x, center_y = result["cell_center"]
            record.update(
                {
                    "picked": {
                        "x": float(center_x),
                        "y": float(center_y),
                        "crs": QgsRasterLayer(
                            self.processor.filled_dem_path, "Filled DEM"
                        ).crs().authid(),
                    },
                    "catchment_area_m2": result["catchment_area_m2"],
                    "mask_path": result["raster_path"],
                    "vector_path": result["vector_path"],
                }
            )

        if state.get("dem_domain"):
            self.prepare_dem_domain(state)
        self.merge_active_domains(state)
        self.update_gauge_domain_ids(state, snapped_layer)
        self.write_gauges(prepared)
        save_state(self.project_folder, state)
        self.remove_deselected_gauges(previously_gauged, state)
        return state

    def validate_gauge_assignments(
        self,
        assignments: Iterable[OutletAssignment],
    ) -> dict[str, _PreparedGauge]:
        """Validate every selected gauge before any workflow output is written."""
        prepared: dict[str, _PreparedGauge] = {}
        numeric_ids: dict[int, str] = {}
        output_folder = Path(streamflow_observation_folder(self.project_folder))
        for assignment in assignments:
            if not assignment.is_gauge:
                continue
            gauge_id = station_id_int(assignment.outlet_id)
            previous = numeric_ids.get(gauge_id)
            if previous is not None:
                raise StationIdError(
                    f"Outlet IDs '{previous}' and '{assignment.outlet_id}' "
                    f"map to the same numeric gauge ID {gauge_id}."
                )
            numeric_ids[gauge_id] = assignment.outlet_id
            layer = assignment.discharge_layer
            if layer is None:
                raise ValueError(
                    f"Select a discharge table for outlet {assignment.outlet_id}."
                )
            is_valid = getattr(layer, "isValid", None)
            if callable(is_valid) and not is_valid():
                raise ValueError(
                    "The discharge table for outlet "
                    f"{assignment.outlet_id} is invalid."
                )
            input_path = local_source_path(layer)
            if input_path is None or input_path.suffix.lower() not in {
                ".csv",
                ".txt",
            }:
                raise ValueError(
                    "The discharge table for outlet "
                    f"{assignment.outlet_id} must be a local CSV or TXT file."
                )
            if not records_from_layer(layer):
                raise ValueError(
                    "No valid discharge records were found for outlet "
                    f"{assignment.outlet_id}."
                )
            output_path = output_folder / streamflow_filename(
                assignment.outlet_id
            )
            prepared[assignment.outlet_id] = _PreparedGauge(
                assignment=assignment,
                gauge_id=gauge_id,
                input_path=input_path.resolve(),
                output_path=output_path,
            )
        return prepared

    def update_gauge_domain_ids(self, state: dict, point_layer=None) -> None:
        """Store IDs of every active domain intersecting each gauge point."""
        assign_domain_ids(state)
        point_layer = point_layer or self.membership_point_layer(state)
        features = self._features_by_outlet(point_layer, require_exact=False)
        domain_geometries = []
        for domain in active_domain_records(state):
            if domain.get("is_dem_domain"):
                path = self.dem_paths()[1]
            else:
                value = domain.get("vector_path")
                path = (
                    str(resolve_output_path(self.project_folder, value))
                    if value
                    else ""
                )
            layer = QgsVectorLayer(path, f"Domain {domain['domain_id']}", "ogr")
            if not path or not layer.isValid():
                raise ValueError(
                    f"Domain polygon is missing for outlet {domain['outlet_id']}."
                )
            geometry = None
            for feature in layer.getFeatures():
                current = QgsGeometry(feature.geometry())
                if current.isEmpty():
                    continue
                geometry = (
                    current if geometry is None else geometry.combine(current)
                )
            if geometry is None or geometry.isEmpty():
                raise ValueError(
                    f"Domain polygon is empty for outlet {domain['outlet_id']}."
                )
            domain_geometries.append(
                (int(domain["domain_id"]), layer.crs(), geometry)
            )

        for outlet_id in gauged_outlet_ids(state):
            feature = features.get(outlet_id)
            if feature is None or feature.geometry().isEmpty():
                raise ValueError(f"Gauge outlet {outlet_id} has no point geometry.")
            memberships = []
            for domain_id, target_crs, domain_geometry in domain_geometries:
                point_geometry = QgsGeometry(feature.geometry())
                source_crs = point_layer.crs()
                if (
                    source_crs.isValid()
                    and target_crs.isValid()
                    and source_crs != target_crs
                ):
                    transform = QgsCoordinateTransform(
                        source_crs,
                        target_crs,
                        QgsProject.instance(),
                    )
                    transform.setBallparkTransformsAreAppropriate(True)
                    point_geometry.transform(transform)
                if domain_geometry.intersects(point_geometry):
                    memberships.append(domain_id)
            if not memberships:
                raise ValueError(
                    f"Gauge outlet {outlet_id} is not inside an active domain."
                )
            state["outlets"][outlet_id]["domain_ids"] = memberships

    def membership_point_layer(self, state: dict):
        """Return snapped gauge points when the saved mode requires them."""
        if state.get("definition_mode") == DOMAIN_MODE_SNAPPED:
            path = getattr(self.processor, "snapped_points_path", None)
            if not path:
                path = os.path.join(
                    geometry_folder(self.project_folder),
                    "2_pour_points_snapped.shp",
                )
            layer = QgsVectorLayer(path, "Snapped pour points", "ogr")
            if layer.isValid():
                return layer
        return self.pour_points_layer

    def regenerate_snapped_points(self) -> str:
        """Recreate snapped points so moved inputs cannot reuse stale positions."""
        if not self.processor._ensure_channel_network(
            self.processor.process_channel_network,
            self.processor.process_flow_accumulation,
            self.processor.fill_dem,
        ):
            raise RuntimeError("The channel network could not be prepared.")
        path = getattr(self.processor, "snapped_points_path", None)
        if not path:
            path = os.path.join(
                geometry_folder(self.project_folder),
                "2_pour_points_snapped.shp",
            )
        if os.path.exists(path):
            self.processor._remove_stale_vector_output(path)
        self.processor.snapped_points_path = None
        if not self.processor._ensure_snapped_points(
            self.processor.snap_points,
            self.processor.process_channel_network,
            self.processor.process_flow_accumulation,
            self.processor.fill_dem,
        ):
            raise RuntimeError("The pour points could not be snapped.")
        return str(self.processor.snapped_points_path)

    @staticmethod
    def validate_unique_state_gauge_ids(state: dict) -> None:
        """Require distinct numeric gauge IDs across all selected gauges."""
        seen = {}
        for outlet_id in gauged_outlet_ids(state):
            gauge_id = station_id_int(outlet_id)
            previous = seen.get(gauge_id)
            if previous is not None:
                raise StationIdError(
                    f"Outlet IDs '{previous}' and '{outlet_id}' map to the "
                    f"same numeric gauge ID {gauge_id}."
                )
            seen[gauge_id] = outlet_id
            state["outlets"][outlet_id]["gauge_id"] = gauge_id

    def prepare_dem_domain(self, state: dict) -> None:
        """Write the domain containing every valid filled-DEM cell."""
        assign_domain_ids(state)
        domain_id = state.get("dem_domain_id")
        if domain_id is None:
            return
        output_raster, output_vector = self.dem_paths()
        reference = self.processor._read_raster_array(
            self.processor.filled_dem_path,
            as_float=True,
        )
        if not reference:
            raise RuntimeError("Could not read the filled DEM.")
        deps = self.processor._get_python_morphology_deps()
        dem, invalid_mask, _ = self.processor._normalise_dem_array(
            reference["array"],
            reference["nodata"],
        )
        if dem is None:
            raise RuntimeError("Could not determine the valid DEM domain.")
        domain = deps["np"].where(~invalid_mask, int(domain_id), 0).astype(
            deps["np"].int32
        )
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

    @staticmethod
    def require_active_domain(state: dict) -> None:
        """Reject a configured workflow without any selected domain."""
        if not active_domain_records(state):
            raise ValueError("Select at least one domain outlet.")

    def merge_active_domains(self, state: dict) -> str | None:
        """Merge every active domain polygon and return the merged path."""
        vector_paths = []
        for domain in active_domain_records(state):
            if domain.get("is_dem_domain"):
                path = self.dem_paths()[1]
            else:
                value = domain.get("vector_path")
                path = (
                    str(resolve_output_path(self.project_folder, value))
                    if value
                    else ""
                )
            if not path or not os.path.exists(path):
                raise ValueError(
                    f"Delineate and save domain {domain['outlet_id']} first."
                )
            layer = QgsVectorLayer(path, Path(path).stem, "ogr")
            if not layer.isValid() or layer.featureCount() < 1:
                raise ValueError(
                    f"Domain polygon is invalid or empty: {Path(path).name}"
                )
            vector_paths.append(path)

        merged_path = os.path.join(
            geometry_folder(self.project_folder),
            "Watersheds",
            "4_watershed_merged_vector.shp",
        )
        if not vector_paths:
            self.processor._remove_vector_dataset(merged_path)
            self.processor.merged_watershed_path = None
            self.processor.watershed_vector_path = None
            return None

        pending_path = os.path.splitext(merged_path)[0] + "_pending.shp"
        self.processor._remove_vector_dataset(pending_path)
        try:
            result = self.processor.run_processing_algorithm(
                "native:mergevectorlayers",
                {"LAYERS": vector_paths, "CRS": None, "OUTPUT": pending_path},
            )
            pending_layer = QgsVectorLayer(
                pending_path,
                "Pending merged domains",
                "ogr",
            )
            valid = (
                bool(result)
                and pending_layer.isValid()
                and pending_layer.featureCount() > 0
            )
            pending_layer = None
            if not valid:
                raise RuntimeError("Could not merge the active domain polygons.")
            self._publish_vector_dataset(pending_path, merged_path)
        finally:
            self.processor._remove_vector_dataset(pending_path)

        self.processor.merged_watershed_path = merged_path
        self.processor.watershed_vector_path = merged_path
        self.processor.mark_output_prepared(
            merged_path,
            name="4_watershed_merged",
            loaded=False,
        )
        return merged_path

    def domain_output_folder(self) -> str:
        path = os.path.join(
            geometry_folder(self.project_folder),
            "Watersheds",
            "DomainDelineation",
        )
        os.makedirs(path, exist_ok=True)
        return path

    def outlet_paths(
        self,
        outlet_id: str,
        state: dict,
        *,
        preview: bool = False,
    ) -> tuple[str, str]:
        prefix = "_preview_" if preview else "4_watershed_"
        ordered = list(state.get("outlet_order", self.outlet_ids))
        try:
            index = ordered.index(str(outlet_id)) + 1
        except ValueError:
            index = len(ordered) + 1
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", str(outlet_id)).strip("_")
        base = os.path.join(
            self.domain_output_folder(),
            f"{prefix}{index}_{safe or 'outlet'}",
        )
        return base + ".tif", base + ".shp"

    def dem_paths(self) -> tuple[str, str]:
        base = os.path.join(self.domain_output_folder(), "4_watershed_DEM")
        return base + ".tif", base + ".shp"

    @staticmethod
    def layer_source(layer) -> str:
        source = str(layer.source() or "")
        local_path = source.split("|", 1)[0]
        if os.path.exists(local_path):
            return str(Path(local_path).resolve())
        return source

    def remove_observation_file(self, outlet_id: str) -> None:
        try:
            filename = streamflow_filename(outlet_id)
        except ValueError as error:
            self.main_dialog.log_message(
                f"WARNING: Could not remove old discharge output: {error}"
            )
            return
        path = Path(streamflow_observation_folder(self.project_folder)) / filename
        try:
            path.unlink(missing_ok=True)
            key = self.processor.output_state_key(str(path))
            self.processor.processing_state.get("outputs", {}).pop(key, None)
            self.processor.save_processing_state()
        except OSError as error:
            self.main_dialog.log_message(
                f"WARNING: Could not remove old discharge output: {error}"
            )

    def _validated_assignment_list(
        self,
        assignments: Iterable[OutletAssignment],
    ) -> list[OutletAssignment]:
        assignments = list(assignments)
        ids = [record.outlet_id for record in assignments]
        if len(ids) != len(set(ids)):
            raise ValueError("Outlet assignments contain duplicate IDs.")
        if ids != self.outlet_ids:
            raise ValueError(
                "Outlet assignments no longer match the selected pour-point layer."
            )
        return assignments

    def apply_assignment_records(
        self,
        state: dict,
        assignments: Iterable[OutletAssignment],
        prepared: dict[str, _PreparedGauge],
    ) -> None:
        for assignment in assignments:
            record = state["outlets"][assignment.outlet_id]
            record["is_gauged"] = bool(assignment.is_gauge)
            record["is_domain"] = bool(assignment.is_domain)
            if not assignment.is_gauge:
                record["discharge_file"] = ""
                record.pop("gauge_id", None)
                record.pop("gauge_filename", None)
                record.pop("gauge_path", None)
                record["domain_ids"] = []
                continue
            gauge = prepared[assignment.outlet_id]
            record.update(
                {
                    "discharge_file": str(gauge.input_path),
                    "gauge_id": gauge.gauge_id,
                    "gauge_filename": gauge.output_path.name,
                    "gauge_path": str(gauge.output_path),
                }
            )

    def write_gauges(self, prepared: dict[str, _PreparedGauge]) -> None:
        for gauge in prepared.values():
            output = write_streamflow_observation(
                gauge.assignment.discharge_layer,
                gauge.assignment.outlet_id,
                gauge.output_path.parent,
            )
            self.processor.mark_output_prepared(
                str(output),
                name=output.name,
                loaded=False,
            )

    def remove_deselected_gauges(self, previous: set[str], state: dict) -> None:
        current = set(gauged_outlet_ids(state))
        for outlet_id in previous - current:
            self.remove_observation_file(outlet_id)

    def _features_by_outlet(
        self,
        layer,
        *,
        require_exact: bool,
    ) -> dict[str, object]:
        field = find_outlet_id_field(layer, self.outlet_id_field)
        features = {}
        for feature in layer.getFeatures():
            outlet_id = station_id_text(feature.attribute(field))
            if not outlet_id:
                raise StationIdError(
                    f"A pour point has an empty {self.outlet_id_field} value."
                )
            if outlet_id in features:
                raise StationIdError(
                    f"Duplicate outlet ID '{outlet_id}' found in snapped pour points."
                )
            features[outlet_id] = feature
        missing = set(self.outlet_ids) - set(features)
        extra = set(features) - set(self.outlet_ids)
        if missing or (require_exact and extra):
            details = []
            if missing:
                details.append(f"missing: {', '.join(sorted(missing))}")
            if require_exact and extra:
                details.append(f"unexpected: {', '.join(sorted(extra))}")
            raise StationIdError(
                "Snapped pour-point IDs do not match the current layer ("
                + "; ".join(details)
                + "). Recreate snapped points."
            )
        return features

    def _polygonize_mask(self, raster_path: str, vector_path: str) -> bool:
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
                and self.processor._copy_nonzero_polygons(raw_path, vector_path)
            )
        finally:
            self.processor._remove_vector_dataset(raw_path)

    def _publish_vector_dataset(self, source_path: str, target_path: str) -> None:
        source_base = os.path.splitext(source_path)[0]
        target_base = os.path.splitext(target_path)[0]
        self.processor._remove_vector_dataset(target_path)
        for extension in (
            ".shp",
            ".shx",
            ".dbf",
            ".prj",
            ".cpg",
            ".qpj",
            ".fix",
        ):
            source = source_base + extension
            if os.path.exists(source):
                os.replace(source, target_base + extension)


__all__ = ["DomainWorkflow"]

# -*- coding: utf-8 -*-
"""Shared domain and gauge assignment workflow."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsGeometry,
    QgsPointXY,
)
from ....core.handlers.raster.tasks import _read_raster, merge_domain_masks
from ....core.handlers.store.paths import domain_data_folder, geometry_folder
from ....core.handlers.store.layout import MERGED_MASK_NAME, domain_dem_path
from ....qt.dialogs.discharge_assignment import OutletAssignment
from ....core.handlers.file.discharge import (
    local_source_path,
    records_from_layer,
    streamflow_filename,
    write_streamflow_observation,
)
from ....core.morphology.hydrology.observation_paths import streamflow_observation_folder
from ....core.morphology.hydrology.outlets import (
    StationIdError,
    find_outlet_id_field,
    station_id_int,
    station_id_text,
)
from ....core.handlers.state.domain_state import (
    DOMAIN_MODE_DEM_EXTENT,
    active_domain_records,
    assign_domain_ids,
    gauged_outlet_ids,
    load_state,
    resolve_output_path,
    save_state,
)
from ... import layers


class _DomainMask:
    """Point-in-domain tests straight off a delineation mask raster."""

    def __init__(self, path: str) -> None:
        raster = _read_raster(path)
        self._values = raster["array"]
        self._transform = raster["geotransform"]
        self.crs = QgsCoordinateReferenceSystem(raster.get("projection") or "")

    def any(self) -> bool:
        return bool((self._values != 0).any())

    def contains(self, point) -> bool:
        origin_x, pixel_width, _rx, origin_y, _ry, pixel_height = self._transform
        col = int((point.x() - origin_x) / pixel_width)
        row = int((point.y() - origin_y) / pixel_height)
        rows, cols = self._values.shape
        if not (0 <= row < rows and 0 <= col < cols):
            return False
        return bool(self._values[row, col] != 0)


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
        point_layer = point_layer or self.pour_points_layer
        missing_points = [
            outlet_id
            for outlet_id in gauged_outlet_ids(state)
            if not isinstance(
                state["outlets"].get(outlet_id, {}).get("gauge_point"), dict
            )
        ]
        features = (
            self._features_by_outlet(point_layer, require_exact=False)
            if missing_points
            else {}
        )
        domain_masks = []
        for domain in active_domain_records(state):
            path = self.domain_mask_path(domain)
            if not path or not os.path.exists(path):
                raise ValueError(
                    f"Domain mask is missing for outlet {domain['outlet_id']}."
                )
            mask = _DomainMask(path)
            if not mask.any():
                raise ValueError(
                    f"Domain mask is empty for outlet {domain['outlet_id']}."
                )
            domain_masks.append((int(domain["domain_id"]), mask))

        for outlet_id in gauged_outlet_ids(state):
            point_geometry, source_crs = self._effective_gauge_geometry(
                state["outlets"][outlet_id],
                features.get(outlet_id),
                point_layer,
            )
            memberships = []
            for domain_id, mask in domain_masks:
                candidate = QgsGeometry(point_geometry)
                transform = layers.transform_between(source_crs, mask.crs)
                if transform is not None:
                    candidate.transform(transform)
                if mask.contains(candidate.asPoint()):
                    memberships.append(domain_id)
            if not memberships:
                raise ValueError(
                    f"Gauge outlet {outlet_id} is not inside an active domain."
                )
            state["outlets"][outlet_id]["domain_ids"] = memberships

    def domain_mask_path(self, domain: dict) -> str:
        """Return the delineation mask raster for one active domain record."""
        if domain.get("is_dem_domain"):
            return self.dem_mask_path()
        value = domain.get("mask_path")
        if not value:
            return ""
        return str(resolve_output_path(self.project_folder, value))

    @staticmethod
    def _effective_gauge_geometry(record, feature, point_layer):
        point = record.get("gauge_point") if isinstance(record, dict) else None
        if isinstance(point, dict):
            try:
                geometry = QgsGeometry.fromPointXY(
                    QgsPointXY(float(point["x"]), float(point["y"]))
                )
            except (KeyError, TypeError, ValueError):
                geometry = None
            if geometry is not None and not geometry.isEmpty():
                source_crs = QgsCoordinateReferenceSystem(
                    str(point.get("crs", "") or "")
                )
                return geometry, source_crs
        if feature is None or feature.geometry().isEmpty():
            raise ValueError("A gauged outlet has no point geometry.")
        return QgsGeometry(feature.geometry()), point_layer.crs()

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
        output_raster = self.dem_mask_path()
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
        state["dem_domain_directory"] = domain_data_folder(
            self.project_folder, "dem_extent"
        )
        state["dem_domain_path"] = domain_dem_path(
            self.project_folder, "dem_extent"
        )

    @staticmethod
    def require_active_domain(state: dict) -> None:
        """Reject a configured workflow without any selected domain."""
        if not active_domain_records(state):
            raise ValueError("Select at least one domain outlet.")

    def merge_active_domains(self, state: dict) -> str | None:
        """Merge every active domain mask and return the merged raster path."""
        masks = []
        for domain in active_domain_records(state):
            path = self.domain_mask_path(domain)
            if not path or not os.path.exists(path):
                raise ValueError(
                    f"Delineate and save domain {domain['outlet_id']} first."
                )
            masks.append((int(domain["domain_id"]), path))

        merged_path = self.merged_mask_path(self.project_folder)
        if not masks:
            self.processor._remove_stale_raster_output(merged_path)
            self.processor.merged_watershed_path = None
            self.processor.watershed_vector_path = None
            return None

        merge_domain_masks(masks, merged_path)
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

    def outlet_mask_path(
        self,
        outlet_id: str,
        state: dict,
        *,
        preview: bool = False,
    ) -> str:
        prefix = "_preview_" if preview else "4_watershed_"
        ordered = list(state.get("outlet_order", self.outlet_ids))
        try:
            index = ordered.index(str(outlet_id)) + 1
        except ValueError:
            index = len(ordered) + 1
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", str(outlet_id)).strip("_")
        return os.path.join(
            self.domain_output_folder(),
            f"{prefix}{index}_{safe or 'outlet'}.tif",
        )

    def dem_mask_path(self) -> str:
        return os.path.join(self.domain_output_folder(), "4_watershed_DEM.tif")

    @staticmethod
    def merged_mask_path(project_folder) -> str:
        return os.path.join(
            geometry_folder(project_folder),
            "Watersheds",
            MERGED_MASK_NAME,
        )

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


__all__ = ["DomainWorkflow"]

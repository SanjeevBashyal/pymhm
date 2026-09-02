# -*- coding: utf-8 -*-
"""Per-pour-point watershed delineation and watershed merge outputs."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    QMessageBox,
    QgsVectorLayer,
    QgsRasterLayer,
    QgsCoordinateTransform,
    QgsProject,
    NULL,
)
from .pour_point_workflow import PourPointWorkflowMixin
from ...qgis_bridge import layers


class WatershedDelineationMixin(PourPointWorkflowMixin):
    """Per-pour-point watershed delineation and watershed merge outputs."""

    def delineate_watershed(self) -> None:
        """Step 4: Delineate upstream basins for snapped points with pyflwdir."""
        self.log_message("\n--- Starting Geometry Step 4: Delineate Watershed ---")
        if not self.check_prerequisites(needs_pour_points=True):
            return

        self.log_message(
            "Preparing snapped pour points before watershed delineation..."
        )
        if not self._ensure_snapped_points(
            self.snap_points,
            self.process_channel_network,
            self.process_flow_accumulation,
            self.fill_dem,
        ):
            return

        if not self._ensure_flow_direction(self.process_flow_direction, self.fill_dem):
            return

        context = self._build_flwdir_from_filled_dem()
        if not context:
            self.log_message("Watershed delineation failed.")
            return

        self.log_message(
            "Delineating watersheds for each snapped point with pyflwdir..."
        )

        snapped_points_layer = QgsVectorLayer(
            self.snapped_points_path, "Snapped Points", "ogr"
        )
        if not snapped_points_layer.isValid():
            QMessageBox.warning(
                self.dialog, "Error", "Could not load snapped points layer."
            )
            return

        features = list(snapped_points_layer.getFeatures())
        if not features:
            QMessageBox.warning(self.dialog, "Error", "No snapped points found.")
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        watershed_output_folder = os.path.join(geometry_folder, "Watersheds")
        os.makedirs(watershed_output_folder, exist_ok=True)
        self.log_message(f"Watershed outputs folder: {watershed_output_folder}")

        watershed_outputs = []
        field_names = snapped_points_layer.fields().names()
        snap_status_field = self._snap_status_field(field_names)
        point_transform = layers.transform_to_raster(
            snapped_points_layer, self.filled_dem_path, log=self.log_message)

        for i, feature in enumerate(features):
            geom = feature.geometry()
            if geom.isEmpty():
                self.log_message(
                    f"Warning: Skipping empty snapped point feature {feature.id()}."
                )
                continue

            point = geom.asPoint()
            name_attr = feature.attribute("Name") if "Name" in field_names else None
            if not name_attr or name_attr == NULL:
                name_attr = f"Watershed_{i + 1}"

            clean_name = "".join(
                c for c in str(name_attr) if c.isalnum() or c in (" ", "-", "_")
            ).rstrip()
            clean_name = clean_name.replace(" ", "_") or f"Watershed_{i + 1}"

            self.log_message(f"Processing watershed for point: {name_attr}")

            if snap_status_field is not None:
                snap_status = feature.attribute(snap_status_field)
                if snap_status in (None, NULL, ""):
                    self.log_message(
                        f"Warning: Skipping point {name_attr}; snapped point has no snap status.")
                    continue
                if str(snap_status) == "failed":
                    self.log_message(
                        f"Warning: Skipping point {name_attr}; point could not be snapped to the channel network.")
                    continue

            point_for_flwdir = point
            if point_transform is not None:
                try:
                    point_for_flwdir = point_transform.transform(point_for_flwdir)
                except Exception as e:
                    self.log_message(
                        f"Warning: Skipping point {name_attr}; could not transform to filled DEM CRS: {e}")
                    continue

            watershed_raster_path = os.path.join(
                watershed_output_folder, f"4_watershed_{clean_name}.tif"
            )
            watershed_vector_path = os.path.join(
                watershed_output_folder, f"4_watershed_{clean_name}.shp"
            )

            try:
                watershed = self.delineate_single_outlet(
                    point_for_flwdir.x(),
                    point_for_flwdir.y(),
                    watershed_raster_path,
                    watershed_vector_path,
                    basin_id=i + 1,
                    context=context,
                )
            except ValueError as e:
                self.log_message(
                    f"Warning: Skipping point {name_attr}; {e}")
                continue

            if not watershed:
                self.log_message(
                    f"Failed to delineate watershed for point: {name_attr}"
                )
                continue

            watershed_outputs.append(
                {
                    "name": name_attr,
                    "clean_name": clean_name,
                    "raster_path": watershed["raster_path"],
                    "vector_path": watershed["vector_path"],
                    "point": point,
                    "catchment_area_m2": watershed["catchment_area_m2"],
                }
            )
            self.load_layer(
                watershed["raster_path"], f"4_Watershed_{clean_name}"
            )

        if not watershed_outputs:
            QMessageBox.warning(
                self.dialog, "Error", "No watersheds were successfully created."
            )
            return

        self.watershed_raster_path = watershed_outputs[0]["raster_path"]
        self.watershed_vector_path = watershed_outputs[0]["vector_path"]

        self.log_message("Merging all watershed vector layers...")
        self.merged_watershed_path = os.path.join(
            watershed_output_folder, "4_watershed_merged_vector.shp"
        )

        vector_layer_paths = [
            watershed_info["vector_path"]
            for watershed_info in watershed_outputs
            if watershed_info.get("vector_path")
            and os.path.exists(watershed_info["vector_path"])
        ]

        if vector_layer_paths:
            self.log_message(
                f"Merging {len(vector_layer_paths)} watershed vector layers..."
            )

            self._remove_vector_dataset(self.merged_watershed_path)
            params_merge = {
                "LAYERS": vector_layer_paths,
                "CRS": None,
                "OUTPUT": self.merged_watershed_path,
            }

            result_merge = self.run_processing_algorithm(
                "native:mergevectorlayers", params_merge
            )

            if result_merge and os.path.exists(self.merged_watershed_path):
                self.watershed_vector_path = self.merged_watershed_path
                self.log_message(
                    f"Merged watershed vector saved: {self.merged_watershed_path}"
                )
            else:
                self.log_message("Warning: Failed to merge watershed vector layers.")
        else:
            self.log_message("Warning: No valid vector layers found to merge.")

        if self.merged_watershed_path and os.path.exists(self.merged_watershed_path):
            self.mark_output_prepared(
                self.merged_watershed_path, name="4_watershed_merged", loaded=False
            )
            self.load_layer(
                self.merged_watershed_path, "4_watershed_merged", is_raster=False
            )
        else:
            self.merged_watershed_path = None

    def delineate_single_outlet(
            self,
            x: float,
            y: float,
            raster_path: str,
            vector_path: str | None = None,
            *,
            basin_id: int = 1,
            context: dict | None = None) -> dict | None:
        """Delineate one outlet whose coordinates use the filled DEM CRS."""
        if not raster_path:
            raise ValueError("A watershed raster output path is required.")

        if context is None:
            if not self._ensure_filled_dem(self.fill_dem):
                return None
            context = self._build_flwdir_from_filled_dem()
            if not context:
                return None

        deps = context["deps"]
        np = deps["np"]
        gdal = deps["gdal"]
        flwdir = context["flwdir"]
        reference = context["reference"]
        try:
            cell = self._point_to_flwdir_cell(x, y, reference, deps)
        except (TypeError, ValueError) as e:
            raise ValueError("Outlet coordinates must be numeric.") from e
        if cell is None:
            raise ValueError(
                f"Outlet coordinates ({float(x):.6f}, {float(y):.6f}) "
                "are outside the filled DEM grid.")

        row, col, center_x, center_y = cell
        if context["invalid_mask"][row, col]:
            raise ValueError(
                f"Outlet coordinates ({float(x):.6f}, {float(y):.6f}) "
                "are outside the valid filled DEM domain.")

        try:
            basin_id = int(basin_id)
        except (TypeError, ValueError) as e:
            raise ValueError("Basin ID must be a positive integer.") from e
        if basin_id < 1:
            raise ValueError("Basin ID must be a positive integer.")

        try:
            basin_map = flwdir.basins(
                xy=([center_x], [center_y]), ids=[basin_id])
        except Exception as e:
            self.log_message(f"Failed to delineate watershed: {e}")
            return None

        basin_mask = basin_map == basin_id
        if not np.any(basin_mask):
            self.log_message("Failed to delineate watershed: basin is empty.")
            return None

        watershed_raster = np.where(
            basin_mask, basin_id, 0).astype(np.int32)
        raster_path = os.path.abspath(raster_path)
        if not self._write_raster_array(
                raster_path,
                watershed_raster,
                reference,
                nodata=0,
                gdal_type=gdal.GDT_Int32):
            return None

        saved_vector_path = None
        if vector_path:
            vector_path = os.path.abspath(vector_path)
            os.makedirs(os.path.dirname(vector_path), exist_ok=True)
            self._remove_vector_dataset(vector_path)
            raw_vector_path = (
                os.path.splitext(vector_path)[0] + "_raw.shp")
            self._remove_vector_dataset(raw_vector_path)
            try:
                result = self.run_processing_algorithm(
                    "gdal:polygonize",
                    {
                        "INPUT": raster_path,
                        "BAND": 1,
                        "FIELD": "DN",
                        "EIGHT_CONNECTEDNESS": False,
                        "EXTRA": "",
                        "OUTPUT": raw_vector_path,
                    },
                )
                if (
                        not result
                        or not os.path.exists(raw_vector_path)
                        or not self._copy_nonzero_polygons(
                            raw_vector_path, vector_path)):
                    self.log_message(
                        "Failed to create the watershed polygon.")
                    return None
                saved_vector_path = vector_path
                self.log_message(
                    f"Watershed polygon saved: {saved_vector_path}")
            finally:
                self._remove_vector_dataset(raw_vector_path)

        area_map = context.get("_upstream_area_m2")
        if area_map is None:
            area_map = flwdir.upstream_area(unit="m2")
            context["_upstream_area_m2"] = area_map
        catchment_area_m2 = float(area_map[row, col])
        if not np.isfinite(catchment_area_m2) or catchment_area_m2 <= 0:
            catchment_area_m2 = (
                int(np.count_nonzero(basin_mask))
                * self._reference_cell_area_m2(reference, deps)
            )

        return {
            "raster_path": raster_path,
            "vector_path": saved_vector_path,
            "cell_center": (center_x, center_y),
            "catchment_area_m2": catchment_area_m2,
        }

    def _snap_status_field(self, field_names: list[str]) -> str | None:
        """Return the snap status field, accounting for shapefile truncation."""
        for candidate in ("snap_status", "snap_statu"):
            if candidate in field_names:
                return candidate
        return None

    def _point_to_flwdir_cell_center(
            self,
            x: float,
            y: float,
            reference: dict,
            deps: dict) -> tuple[float, float] | None:
        """Return the raster cell center for a map coordinate, or None if outside."""
        cell = self._point_to_flwdir_cell(x, y, reference, deps)
        if cell is None:
            return None
        return cell[2], cell[3]

    def _point_to_flwdir_cell(
            self,
            x: float,
            y: float,
            reference: dict,
            deps: dict) -> tuple[int, int, float, float] | None:
        """Return row, column, and cell center for a map coordinate."""
        transform = deps["Affine"].from_gdal(*reference["geotransform"])
        col_float, row_float = (~transform) * (float(x), float(y))
        np = deps["np"]
        row = int(np.floor(row_float))
        col = int(np.floor(col_float))
        if row < 0 or col < 0:
            return None
        if row >= int(reference["rows"]) or col >= int(reference["cols"]):
            return None

        center_x, center_y = transform * (col + 0.5, row + 0.5)
        return row, col, float(center_x), float(center_y)

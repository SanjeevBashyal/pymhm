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

        deps = context["deps"]
        np = deps["np"]
        gdal = deps["gdal"]
        flwdir = context["flwdir"]
        reference = context["reference"]

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
        point_transform = self._snapped_to_filled_dem_transform(snapped_points_layer)

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

            flwdir_xy = self._point_to_flwdir_cell_center(
                point_for_flwdir.x(),
                point_for_flwdir.y(),
                reference,
                deps,
            )
            if flwdir_xy is None:
                self.log_message(
                    f"Warning: Skipping point {name_attr}; snapped coordinates "
                    f"({point_for_flwdir.x():.6f}, {point_for_flwdir.y():.6f}) "
                    "are outside the filled DEM domain.")
                continue

            watershed_raster_path = os.path.join(
                watershed_output_folder, f"4_watershed_{clean_name}.tif"
            )
            watershed_vector_raw_path = os.path.join(
                watershed_output_folder, f"4_watershed_{clean_name}_raw.shp"
            )
            watershed_vector_path = os.path.join(
                watershed_output_folder, f"4_watershed_{clean_name}.shp"
            )

            basin_id = i + 1
            try:
                basin_map = flwdir.basins(
                    xy=([flwdir_xy[0]], [flwdir_xy[1]]), ids=[basin_id])
            except Exception as e:
                self.log_message(f"Failed to delineate watershed for {name_attr}: {e}")
                continue

            watershed_raster = np.where(basin_map == basin_id, basin_id, 0).astype(
                np.int32
            )
            wrote_raster = self._write_raster_array(
                watershed_raster_path,
                watershed_raster,
                reference,
                nodata=0,
                gdal_type=gdal.GDT_Int32,
            )
            if not wrote_raster:
                self.log_message(
                    f"Failed to write watershed raster for point: {name_attr}"
                )
                continue

            self.log_message(f"Polygonizing watershed raster for: {name_attr}")
            self._remove_vector_dataset(watershed_vector_raw_path)
            params_poly = {
                "INPUT": watershed_raster_path,
                "BAND": 1,
                "FIELD": "DN",
                "EIGHT_CONNECTEDNESS": False,
                "EXTRA": "",
                "OUTPUT": watershed_vector_raw_path,
            }

            result_poly = self.run_processing_algorithm("gdal:polygonize", params_poly)

            if result_poly and os.path.exists(watershed_vector_raw_path):
                if self._copy_nonzero_polygons(
                    watershed_vector_raw_path, watershed_vector_path
                ):
                    self.log_message(
                        f"Watershed polygon saved: {watershed_vector_path}"
                    )
                    watershed_outputs.append(
                        {
                            "name": name_attr,
                            "clean_name": clean_name,
                            "raster_path": watershed_raster_path,
                            "vector_path": watershed_vector_path,
                            "point": point,
                        }
                    )
                    self.load_layer(watershed_raster_path, f"4_Watershed_{clean_name}")
                else:
                    self.log_message(
                        f"Warning: Failed to filter watershed polygon for: {name_attr}"
                    )

                self._remove_vector_dataset(watershed_vector_raw_path)
            else:
                self.log_message(
                    f"Warning: Failed to polygonize watershed for: {name_attr}"
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

    def _snap_status_field(self, field_names: list[str]) -> str | None:
        """Return the snap status field, accounting for shapefile truncation."""
        for candidate in ("snap_status", "snap_statu"):
            if candidate in field_names:
                return candidate
        return None

    def _snapped_to_filled_dem_transform(self, snapped_layer: object) -> object | None:
        """Return a coordinate transform when snapped points and filled DEM CRS differ."""
        filled_dem_layer = QgsRasterLayer(self.filled_dem_path, "Filled_DEM")
        if not filled_dem_layer.isValid():
            return None

        source_crs = snapped_layer.crs()
        target_crs = filled_dem_layer.crs()
        if not source_crs.isValid() or not target_crs.isValid():
            return None
        if source_crs.authid() == target_crs.authid():
            return None

        transform = QgsCoordinateTransform(
            source_crs,
            target_crs,
            QgsProject.instance(),
        )
        transform.setBallparkTransformsAreAppropriate(True)
        self.log_message(
            f"Transforming snapped points from {source_crs.authid()} "
            f"to filled DEM CRS {target_crs.authid()} for watershed delineation.")
        return transform

    def _point_to_flwdir_cell_center(
            self,
            x: float,
            y: float,
            reference: dict,
            deps: dict) -> tuple[float, float] | None:
        """Return the raster cell center for a map coordinate, or None if outside."""
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
        return float(center_x), float(center_y)

# -*- coding: utf-8 -*-
"""Elevation-band raster generation and area summary writing."""
from ..common import (
    os,
    project_geometry_folder,
    csv,
    QMessageBox,
)


class ElevationBandRasterMixin:
    """Elevation-band raster generation and area summary writing."""

    def process_elevation_bands(self):
        """
        Create elevation-band rasters per subcatchment and a CSV area summary.

        The band range is rounded to a clean interval and applied separately to
        each delineated subcatchment raster in Geometry/Watersheds.
        """
        self.log_message("\n--- Creating Elevation Bands ---")
        if not self.check_prerequisites(needs_pour_points=True):
            return

        dem_layer = self.dialog.mMapLayerComboBox_dem.currentLayer()
        if not dem_layer:
            QMessageBox.warning(
                self.dialog, "Input Error", "Please select a DEM Raster Layer.")
            return

        deps = self._get_python_morphology_deps()
        if not deps:
            return

        np = deps["np"]
        gdal = deps["gdal"]

        selected_dem_reference = self._read_raster_array(
            dem_layer.source(), as_float=True)
        if not selected_dem_reference:
            return

        selected_dem_array = selected_dem_reference["array"]
        selected_dem_valid = self._valid_raster_mask(
            selected_dem_array, selected_dem_reference["nodata"], np)
        if not np.any(selected_dem_valid):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "The selected DEM has no valid elevation values.")
            return

        selected_min = float(np.nanmin(selected_dem_array[selected_dem_valid]))
        selected_max = float(np.nanmax(selected_dem_array[selected_dem_valid]))
        band_width = self._ask_elevation_band_width(selected_min, selected_max)
        if band_width is None:
            self.log_message("Elevation band creation cancelled by user.")
            return

        if not self._ensure_filled_dem():
            return
        if not self._ensure_merged_watershed():
            return

        dem_reference = self._read_raster_array(self.filled_dem_path, as_float=True)
        if not dem_reference:
            return

        dem_array = dem_reference["array"]
        dem_valid = self._valid_raster_mask(
            dem_array, dem_reference["nodata"], np)
        if not np.any(dem_valid):
            QMessageBox.warning(
                self.dialog, "Input Error",
                "The filled DEM has no valid elevation values.")
            return

        dem_min = float(np.nanmin(dem_array[dem_valid]))
        dem_max = float(np.nanmax(dem_array[dem_valid]))
        elevation_range = self._elevation_range_from_window_width(
            dem_min, dem_max, band_width)
        if not elevation_range:
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please enter an elevation window width greater than zero.")
            return

        nice_min, nice_max, band_count = elevation_range
        band_edges = np.array(
            [nice_min + i * band_width for i in range(band_count + 1)],
            dtype=np.float64
        )

        self.log_message(
            f"Selected DEM elevation range: {selected_min:.2f} to {selected_max:.2f}")
        self.log_message(
            f"Rounded band range used: {nice_min:.2f} to {nice_max:.2f} "
            f"({band_count} bands, window width {band_width:.2f})")

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        watershed_output_folder = os.path.join(geometry_folder, "Watersheds")
        watershed_rasters = self._collect_watershed_rasters(watershed_output_folder)

        if not watershed_rasters:
            self.log_message(
                "No per-subcatchment watershed rasters found. Running Watershed Delineation first...")
            self.without_layer_loading(self.delineate_watershed)
            watershed_rasters = self._collect_watershed_rasters(watershed_output_folder)

        if not watershed_rasters:
            QMessageBox.warning(
                self.dialog, "Missing Watersheds",
                "No subcatchment watershed rasters were found or created.")
            return

        output_folder = os.path.join(geometry_folder, "ElevationBands")
        os.makedirs(output_folder, exist_ok=True)
        csv_path = os.path.join(output_folder, "elevation_band_areas.csv")
        cell_area_m2 = self._reference_cell_area_m2(dem_reference, deps)

        csv_rows = []
        created_count = 0

        full_band_ids = np.digitize(
            dem_array, band_edges[1:-1], right=False).astype(np.int16) + 1
        full_band_ids = np.clip(full_band_ids, 1, band_count)

        for watershed_path in watershed_rasters:
            watershed_reference = self._read_raster_array(watershed_path)
            if not watershed_reference:
                continue

            watershed_array = watershed_reference["array"]
            if watershed_array.shape != dem_array.shape:
                self.log_message(
                    f"WARNING: Skipping watershed raster with mismatched shape: {watershed_path}")
                continue

            subcatchment_name = os.path.splitext(
                os.path.basename(watershed_path))[0]
            if subcatchment_name.startswith("4_watershed_"):
                subcatchment_name = subcatchment_name[len("4_watershed_"):]
            clean_name = self._clean_output_name(
                subcatchment_name, f"subcatchment_{created_count + 1}")

            watershed_nodata = watershed_reference["nodata"]
            if watershed_nodata is not None and np.isfinite(watershed_nodata):
                watershed_mask = ~np.isclose(watershed_array, watershed_nodata)
            else:
                watershed_mask = np.ones(watershed_array.shape, dtype=bool)
            watershed_mask &= watershed_array > 0

            valid_subcatchment = watershed_mask & dem_valid
            if not np.any(valid_subcatchment):
                self.log_message(
                    f"WARNING: No valid DEM cells inside subcatchment: {clean_name}")
                continue

            elevation_band_array = np.zeros(dem_array.shape, dtype=np.int16)
            elevation_band_array[valid_subcatchment] = full_band_ids[valid_subcatchment]

            output_raster_path = os.path.join(
                output_folder, f"elevation_bands_{clean_name}.tif")
            if self._write_raster_array(
                    output_raster_path,
                    elevation_band_array,
                    dem_reference,
                    nodata=0,
                    gdal_type=gdal.GDT_Int16):
                self.load_layer(
                    output_raster_path,
                    f"Elevation_Bands_{clean_name}")
                created_count += 1
            else:
                self.log_message(
                    f"ERROR: Failed to write elevation bands raster for {clean_name}.")
                continue

            for band_id in range(1, band_count + 1):
                band_mask = valid_subcatchment & (elevation_band_array == band_id)
                area_cells = int(np.count_nonzero(band_mask))
                area_m2 = area_cells * cell_area_m2
                csv_rows.append({
                    "subcatchment": clean_name,
                    "band_id": band_id,
                    "min_elevation": band_edges[band_id - 1],
                    "max_elevation": band_edges[band_id],
                    "area_cells": area_cells,
                    "area_m2": area_m2,
                    "area_km2": area_m2 / 1000000.0
                })

        if created_count == 0:
            QMessageBox.warning(
                self.dialog, "No Output",
                "No elevation band rasters were created.")
            return

        with open(csv_path, "w", newline="", encoding="utf-8") as csv_file:
            fieldnames = [
                "subcatchment",
                "band_id",
                "min_elevation",
                "max_elevation",
                "area_cells",
                "area_m2",
                "area_km2"
            ]
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            for row in csv_rows:
                writer.writerow({
                    "subcatchment": row["subcatchment"],
                    "band_id": row["band_id"],
                    "min_elevation": f"{row['min_elevation']:.6f}",
                    "max_elevation": f"{row['max_elevation']:.6f}",
                    "area_cells": row["area_cells"],
                    "area_m2": f"{row['area_m2']:.6f}",
                    "area_km2": f"{row['area_km2']:.6f}"
                })

        self.log_message(f"Elevation band CSV written: {csv_path}")
        self.mark_output_prepared(
            csv_path,
            name="elevation_band_areas.csv",
            loaded=False
        )
        QMessageBox.information(
            self.dialog,
            "Elevation Bands Created",
            f"Created {created_count} elevation band raster(s).\n{csv_path}")

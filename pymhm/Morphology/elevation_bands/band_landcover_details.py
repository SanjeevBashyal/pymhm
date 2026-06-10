# -*- coding: utf-8 -*-
"""Elevation-band by land-cover detail CSV generation."""
from __future__ import annotations

from ..common import (
    os,
    project_geometry_folder,
    csv,
    QMessageBox,
    QgsRasterLayer,
)
from ..layers.land_cover_class_names import LandCoverClassNameMixin
from .band_landcover_helpers import BandLandCoverHelperMixin
from .band_summary import BandSummaryMixin
from .elevation_band_rasters import ElevationBandRasterMixin
from .raster_alignment import RasterAlignmentMixin
from .raster_masks import RasterMaskMixin


class BandLandCoverDetailsMixin(
        ElevationBandRasterMixin,
        BandSummaryMixin,
        RasterAlignmentMixin,
        BandLandCoverHelperMixin,
        LandCoverClassNameMixin,
        RasterMaskMixin):
    """Elevation-band by land-cover detail CSV generation."""

    def process_band_details(self) -> None:
        """
        Intersect elevation bands with raw land-cover classes and write details.

        The output table has one row per subcatchment and band, plus one dynamic
        area column for each land-cover class found inside the elevation bands.
        Dynamic class-area columns are reported in square metres and use class
        names from the land-cover lookup layer when available.
        """
        self.log_message("\n--- Creating Elevation Band Land Cover Details ---")

        if not self.dialog.project_folder:
            QMessageBox.critical(
                self.dialog,
                "Missing Input",
                "Please select a project folder before proceeding.")
            return

        geometry_folder = project_geometry_folder(self.dialog.project_folder)
        elevation_band_folder = os.path.join(geometry_folder, "ElevationBands")
        elevation_band_summary_path = os.path.join(
            elevation_band_folder, "elevation_band_areas.csv")

        elevation_band_rasters = self._collect_elevation_band_rasters(
            elevation_band_folder)
        if not elevation_band_rasters or not os.path.exists(elevation_band_summary_path):
            self.log_message(
                "Elevation-band rasters or summary CSV are missing. Running Elevation Bands first...")
            self.without_layer_loading(self.process_elevation_bands)
            elevation_band_rasters = self._collect_elevation_band_rasters(
                elevation_band_folder)

        if not elevation_band_rasters:
            QMessageBox.warning(
                self.dialog,
                "Missing Elevation Bands",
                "No elevation band rasters were found or created.")
            return

        land_cover_layer = self.dialog.mMapLayerComboBox_land_cover.currentLayer()
        if not land_cover_layer or not isinstance(land_cover_layer, QgsRasterLayer):
            QMessageBox.warning(
                self.dialog,
                "Input Error",
                "Please select the original Land Cover raster layer for band details.")
            return

        deps = self._get_python_morphology_deps()
        if not deps:
            return

        np = deps["np"]

        first_band_reference = self._read_raster_array(elevation_band_rasters[0])
        if not first_band_reference:
            return

        land_cover_path = self._source_path_for_gdal(land_cover_layer)
        self.log_message(
            f"Using raw land cover classes from: {land_cover_layer.name()} | {land_cover_path}")
        class_names = self._read_land_cover_class_names()
        if class_names:
            self.log_message(
                f"Loaded {len(class_names)} land cover class name(s) for CSV headers.")

        land_use_reference = self._read_raster_array(land_cover_path, as_float=True)
        if not land_use_reference:
            return

        if not self._raster_matches_reference(
                land_use_reference, first_band_reference):
            aligned_land_use_path = os.path.join(
                elevation_band_folder, "land_cover_classes_aligned_to_elevation_bands.tif")
            self.log_message(
                "Land-cover raster grid differs from elevation bands. Aligning with nearest-neighbour resampling...")
            if not self._align_raster_to_reference(
                    land_cover_path,
                    aligned_land_use_path,
                    first_band_reference,
                    deps):
                QMessageBox.critical(
                    self.dialog,
                    "Alignment Error",
                    "Could not align land-use raster to the elevation-band grid.")
                return
            land_use_reference = self._read_raster_array(
                aligned_land_use_path, as_float=True)
            if not land_use_reference:
                return

        land_use_array = land_use_reference["array"]
        land_use_valid = self._valid_raster_mask(
            land_use_array, land_use_reference["nodata"], np)

        summary_rows = self._load_elevation_band_summary(
            elevation_band_summary_path)
        cell_area_m2 = self._reference_cell_area_m2(first_band_reference, deps)

        band_references = []
        class_values = {}
        for band_raster_path in elevation_band_rasters:
            band_reference = self._read_raster_array(band_raster_path)
            if not band_reference:
                continue

            if (band_reference["rows"] != land_use_array.shape[0] or
                    band_reference["cols"] != land_use_array.shape[1]):
                self.log_message(
                    f"WARNING: Skipping band raster with mismatched shape: {band_raster_path}")
                continue

            band_array = band_reference["array"]
            band_valid = self._valid_raster_mask(
                band_array, band_reference["nodata"], np)
            band_valid &= band_array > 0
            valid_intersection = band_valid & land_use_valid

            if np.any(valid_intersection):
                for class_value in np.unique(land_use_array[valid_intersection]):
                    canonical_value = self._canonical_land_cover_class(class_value)
                    class_values[canonical_value] = canonical_value

            subcatchment = self._subcatchment_from_elevation_band_path(
                band_raster_path)
            band_references.append((subcatchment, band_reference))

        if not band_references:
            QMessageBox.warning(
                self.dialog,
                "No Valid Bands",
                "No elevation band rasters matched the land-use raster grid.")
            return

        if not class_values:
            QMessageBox.warning(
                self.dialog,
                "No Land Cover Classes",
                "No valid land-cover classes were found inside the elevation bands.")
            return

        sorted_class_values = sorted(class_values.values(), key=float)
        class_headers = [
            self._land_cover_class_header(class_value, class_names)
            for class_value in sorted_class_values
        ]

        fieldnames = [
            "subcatchment",
            "band_id",
            "min_elevation",
            "max_elevation",
            "area_cells",
            "area_m2",
            "area_km2",
        ] + class_headers

        output_csv_path = os.path.join(
            elevation_band_folder, "elevation_band_land_cover_areas.csv")
        output_rows = []

        for subcatchment, band_reference in band_references:
            band_array = band_reference["array"]
            band_valid = self._valid_raster_mask(
                band_array, band_reference["nodata"], np)
            band_valid &= band_array > 0

            summary_band_ids = sorted(
                band_id for (summary_subcatchment, band_id) in summary_rows
                if summary_subcatchment == subcatchment
            )
            raster_band_ids = sorted(
                int(value) for value in np.unique(band_array[band_valid])
            ) if np.any(band_valid) else []
            band_ids = summary_band_ids or raster_band_ids

            for band_id in band_ids:
                band_mask = band_valid & np.isclose(band_array, band_id)
                summary = summary_rows.get((subcatchment, band_id), {})

                area_cells = summary.get("area_cells")
                if area_cells is None:
                    area_cells = int(np.count_nonzero(band_mask))

                area_m2 = summary.get("area_m2")
                if area_m2 is None:
                    area_m2 = area_cells * cell_area_m2

                area_km2 = summary.get("area_km2")
                if area_km2 is None:
                    area_km2 = area_m2 / 1000000.0

                row = {
                    "subcatchment": subcatchment,
                    "band_id": band_id,
                    "min_elevation": self._format_optional_float(
                        summary.get("min_elevation")),
                    "max_elevation": self._format_optional_float(
                        summary.get("max_elevation")),
                    "area_cells": area_cells,
                    "area_m2": f"{area_m2:.6f}",
                    "area_km2": f"{area_km2:.6f}",
                }

                land_cover_band_mask = band_mask & land_use_valid
                for class_value, class_header in zip(
                        sorted_class_values, class_headers):
                    class_mask = land_cover_band_mask & np.isclose(
                        land_use_array, float(class_value))
                    class_area_m2 = int(np.count_nonzero(class_mask)) * cell_area_m2
                    row[class_header] = f"{class_area_m2:.6f}"

                output_rows.append(row)

        if not output_rows:
            QMessageBox.warning(
                self.dialog,
                "No Output",
                "No elevation-band land-cover detail rows were created.")
            return

        with open(output_csv_path, "w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_rows)

        self.log_message(
            f"Elevation band land-cover details CSV written: {output_csv_path}")
        self.mark_output_prepared(
            output_csv_path,
            name="elevation_band_land_cover_areas.csv",
            loaded=False
        )
        QMessageBox.information(
            self.dialog,
            "Band Details Created",
            f"Created {len(output_rows)} band detail row(s).\n{output_csv_path}")

# -*- coding: utf-8 -*-
"""Elevation-band raster discovery and summary CSV parsing."""
from __future__ import annotations

from ..common import (
    os,
    csv,
)
from ..core.naming import NamingAndRangeMixin


class BandSummaryMixin(NamingAndRangeMixin):
    """Elevation-band raster discovery and summary CSV parsing."""

    def _collect_elevation_band_rasters(
            self,
            elevation_band_folder: str) -> list[str]:
        """Find per-subcatchment elevation-band rasters."""
        if not os.path.exists(elevation_band_folder):
            return []

        raster_paths = []
        for filename in os.listdir(elevation_band_folder):
            lower_name = filename.lower()
            if lower_name.startswith("elevation_bands_") and lower_name.endswith(".tif"):
                raster_paths.append(os.path.join(elevation_band_folder, filename))

        return sorted(raster_paths)

    def _subcatchment_from_elevation_band_path(self, raster_path: str) -> str:
        """Extract the subcatchment name used by elevation-band outputs."""
        stem = os.path.splitext(os.path.basename(raster_path))[0]
        prefix = "elevation_bands_"
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
        return self._clean_output_name(stem, "subcatchment")

    def _load_elevation_band_summary(
            self,
            csv_path: str) -> dict[tuple[str, int], dict[str, float | int | None]]:
        """Load elevation band min/max and area summary rows."""
        if not os.path.exists(csv_path):
            return {}

        summary = {}
        try:
            with open(csv_path, "r", newline="", encoding="utf-8") as csv_file:
                reader = csv.DictReader(csv_file)
                for row in reader:
                    subcatchment = str(row.get("subcatchment", "")).strip()
                    if not subcatchment:
                        continue
                    try:
                        band_id = int(float(row.get("band_id", "")))
                    except (TypeError, ValueError):
                        continue

                    summary[(subcatchment, band_id)] = {
                        "min_elevation": self._parse_optional_float(
                            row.get("min_elevation")),
                        "max_elevation": self._parse_optional_float(
                            row.get("max_elevation")),
                        "area_cells": self._parse_optional_int(
                            row.get("area_cells")),
                        "area_m2": self._parse_optional_float(
                            row.get("area_m2")),
                        "area_km2": self._parse_optional_float(
                            row.get("area_km2")),
                    }
        except Exception as e:
            self.log_message(
                f"WARNING: Could not read elevation band summary CSV: {e}")

        return summary

    def _parse_optional_float(self, value: object) -> float | None:
        """Parse a float value from CSV text, preserving blanks as None."""
        if value in (None, ""):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _parse_optional_int(self, value: object) -> int | None:
        """Parse an integer value from CSV text, preserving blanks as None."""
        if value in (None, ""):
            return None
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return None

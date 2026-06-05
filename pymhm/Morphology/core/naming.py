# -*- coding: utf-8 -*-
"""Output naming, project geometry paths, and elevation range utilities."""
from ..common import (
    os,
    math,
    geometry_folder,
)


class NamingAndRangeMixin:
    """Output naming, project geometry paths, and elevation range utilities."""

    def _clean_output_name(self, value, fallback):
        """Return a filesystem-safe short name."""
        clean_name = "".join(
            c for c in str(value)
            if c.isalnum() or c in (' ', '-', '_')
        ).rstrip()
        return clean_name.replace(' ', '_') or fallback

    def _nice_step(self, raw_step):
        """Round a raw interval up to a 1/2/5/10 style step."""
        if raw_step <= 0:
            return 1.0

        exponent = math.floor(math.log10(raw_step))
        base = 10 ** exponent
        for multiplier in (1, 2, 5, 10):
            step = multiplier * base
            if step >= raw_step:
                return float(step)

        return float(10 * base)

    def _next_nice_step(self, step):
        """Return the next larger 1/2/5/10 style step."""
        if step <= 0:
            return 1.0

        exponent = math.floor(math.log10(step))
        base = 10 ** exponent
        normalized = step / base

        if normalized < 2:
            return float(2 * base)
        if normalized < 5:
            return float(5 * base)
        if normalized < 10:
            return float(10 * base)
        return float(20 * base)

    def _nice_elevation_range(self, min_value, max_value, band_count):
        """Build rounded min/max and edges for a fixed number of elevation bands."""
        if band_count < 1:
            band_count = 1

        if max_value <= min_value:
            max_value = min_value + 1.0

        raw_step = (max_value - min_value) / float(band_count)
        step = self._nice_step(raw_step)
        nice_min = math.floor(min_value / step) * step
        nice_max = nice_min + band_count * step

        while nice_max < max_value:
            step = self._next_nice_step(step)
            nice_min = math.floor(min_value / step) * step
            nice_max = nice_min + band_count * step

        return nice_min, nice_max, step

    def _elevation_range_from_window_width(self, min_value, max_value, window_width):
        """Build rounded elevation range from a requested band window width."""
        if window_width <= 0:
            return None

        if max_value <= min_value:
            max_value = min_value + window_width

        rounded_min = math.floor(min_value / window_width) * window_width
        rounded_max = math.ceil(max_value / window_width) * window_width

        if rounded_max <= rounded_min:
            rounded_max = rounded_min + window_width

        band_count = int(math.ceil(
            ((rounded_max - rounded_min) / window_width) - 1e-12))
        rounded_max = rounded_min + band_count * window_width

        return rounded_min, rounded_max, band_count

    def _geometry_path(self, filename):
        """Return an output path in the project Geometry folder."""
        return os.path.join(geometry_folder(self.dialog.project_folder), filename)

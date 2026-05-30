# -*- coding: utf-8 -*-
"""Reusable raster nodata/validity mask helpers."""


class RasterMaskMixin:
    """Reusable raster nodata/validity mask helpers."""

    def _valid_raster_mask(self, array, nodata, np):
        """Return True where raster values are finite and not nodata."""
        valid_mask = np.isfinite(array)
        if nodata is not None and np.isfinite(nodata):
            valid_mask &= ~np.isclose(array, nodata)
        return valid_mask

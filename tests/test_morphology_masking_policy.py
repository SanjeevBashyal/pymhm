"""Which morphology layers are padded, and which are watershed-masked."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from mhm_qgis import standalone  # noqa: E402

standalone.install(force=True)

from mhm_qgis.core.morphology.layers import CATEGORICAL_PAD_VALUE  # noqa: E402
from mhm_qgis.core.executions.morphology.setup import raster_plan  # noqa: E402


# Layers the DEM workflow produces, and the policy each one must follow.
DEM_DERIVATIVES = (
    "1_dem_filled.tif",
    "1_dem_aspect.tif",
    "1_dem_slope.tif",
    "2_flow_accumulation.tif",
    "2_flow_direction.tif",
    "2_gauge_position.tif",
)
CLASS_LAYERS = (
    "3_land_use.tif",
    "3_soil.tif",
    "3_geology_processed.tif",
)


def _geometry(tmp_path, names):
    folder = tmp_path / "mhm-plugin" / "Z Temp" / "Geometry"
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        (folder / name).touch()
    return str(folder)


def _by_file(layers):
    return {Path(entry["input"]).name: entry for entry in layers}


def test_only_class_layers_pad_with_a_class(tmp_path):
    _geometry(tmp_path, DEM_DERIVATIVES + CLASS_LAYERS)

    layers = _by_file(raster_plan(tmp_path))

    assert set(layers) == set(DEM_DERIVATIVES + CLASS_LAYERS)
    for name in CLASS_LAYERS:
        assert layers[name]["pad"] == CATEGORICAL_PAD_VALUE
    for name in DEM_DERIVATIVES:
        assert layers[name]["pad"] is None


def test_only_dem_derivatives_are_watershed_masked(tmp_path):
    """Class layers keep the full model extent; masking them would cut holes."""
    _geometry(tmp_path, DEM_DERIVATIVES + CLASS_LAYERS)
    layers = {
        Path(entry["masked"]).name: entry for entry in raster_plan(tmp_path)
    }

    assert len(layers) == len(DEM_DERIVATIVES) + len(CLASS_LAYERS)
    for name in CLASS_LAYERS:
        masked = name.replace(".tif", "_masked.tif")
        assert layers[masked]["watershed_masked"] is False
        assert layers[masked]["pad"] == CATEGORICAL_PAD_VALUE
    for name in DEM_DERIVATIVES:
        masked = name.replace(".tif", "_masked.tif")
        assert layers[masked]["watershed_masked"] is True

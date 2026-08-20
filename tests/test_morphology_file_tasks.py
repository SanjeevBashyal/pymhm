from __future__ import annotations

import os

import numpy as np
import pytest
from osgeo import gdal, ogr, osr

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from pymhm.Morphology.file_tasks import (  # noqa: E402
    crop_aligned_l0_raster,
    delineate_domains_file,
    delineate_outlet_file,
    fill_dem_file,
    hydrology_files,
    mask_aligned_l0_raster,
    materialize_domain_dem_file,
    terrain_files,
)


def _dem(path):
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path), 10, 10, 1, gdal.GDT_Float32
    )
    dataset.SetGeoTransform((0, 100, 0, 1000, 0, -100))
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(32632)
    dataset.SetProjection(crs.ExportToWkt())
    values = np.add.outer(
        np.arange(10, 0, -1), np.arange(10, 0, -1)
    ).astype("float32")
    dataset.GetRasterBand(1).WriteArray(values)
    dataset.GetRasterBand(1).SetNoDataValue(-9999)
    dataset = None


def test_aligned_l0_crop_preserves_values_and_pads_with_nodata(tmp_path):
    source = tmp_path / "dem.tif"
    output = tmp_path / "crop.tif"
    _dem(source)

    crop_aligned_l0_raster(
        source,
        output,
        {
            "ncols": 12,
            "nrows": 12,
            "xllcorner": -100.0,
            "yllcorner": -100.0,
            "cellsize": 100.0,
            "nodata_value": -9999.0,
        },
        reference_path=source,
    )

    source_values = gdal.Open(str(source)).ReadAsArray()
    result = gdal.Open(str(output))
    values = result.ReadAsArray()
    assert values.shape == (12, 12)
    np.testing.assert_array_equal(values[1:11, 1:11], source_values)
    assert np.all(values[0] == -9999)
    assert np.all(values[:, 0] == -9999)


def _polygon(path, xmin, xmax, ymin, ymax):
    """Write a single-rectangle shapefile in the DEM CRS."""
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(32632)
    source = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(str(path))
    layer = source.CreateLayer("mask", crs, ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in (
        (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)
    ):
        ring.AddPoint(x, y)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
    source = None


def test_aligned_l0_mask_keeps_values_inside_and_nodata_outside(tmp_path):
    source = tmp_path / "dem.tif"
    mask_vector = tmp_path / "mask.shp"
    output = tmp_path / "masked.tif"
    _dem(source)
    # Cover exactly the four cells spanning columns 2-3 and rows 2-3.
    _polygon(mask_vector, 200.0, 400.0, 600.0, 800.0)

    mask_aligned_l0_raster(
        source,
        output,
        {
            "ncols": 10,
            "nrows": 10,
            "xllcorner": 0.0,
            "yllcorner": 0.0,
            "cellsize": 100.0,
            "nodata_value": -9999.0,
        },
        mask_vector,
        reference_path=source,
    )

    source_values = gdal.Open(str(source)).ReadAsArray()
    values = gdal.Open(str(output)).ReadAsArray()
    assert values.shape == source_values.shape
    np.testing.assert_array_equal(values[2:4, 2:4], source_values[2:4, 2:4])
    kept = np.zeros(values.shape, dtype=bool)
    kept[2:4, 2:4] = True
    assert np.all(values[~kept] == -9999)


def test_aligned_l0_mask_pads_an_expanded_extent_with_nodata(tmp_path):
    source = tmp_path / "dem.tif"
    mask_vector = tmp_path / "mask.shp"
    output = tmp_path / "masked.tif"
    _dem(source)
    _polygon(mask_vector, -200.0, 1200.0, -200.0, 1200.0)

    mask_aligned_l0_raster(
        source,
        output,
        {
            "ncols": 12,
            "nrows": 12,
            "xllcorner": -100.0,
            "yllcorner": -100.0,
            "cellsize": 100.0,
            "nodata_value": -9999.0,
        },
        mask_vector,
        reference_path=source,
    )

    source_values = gdal.Open(str(source)).ReadAsArray()
    values = gdal.Open(str(output)).ReadAsArray()
    assert values.shape == (12, 12)
    np.testing.assert_array_equal(values[1:11, 1:11], source_values)
    assert np.all(values[0] == -9999)
    assert np.all(values[-1] == -9999)


def test_misaligned_l0_source_is_rejected_rather_than_resampled(tmp_path):
    source = tmp_path / "shifted.tif"
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(source), 10, 10, 1, gdal.GDT_Float32
    )
    dataset.SetGeoTransform((50, 100, 0, 1000, 0, -100))
    dataset.GetRasterBand(1).WriteArray(np.ones((10, 10), dtype="float32"))
    dataset = None

    with pytest.raises(ValueError, match="aligned"):
        crop_aligned_l0_raster(
            source,
            tmp_path / "crop.tif",
            {
                "ncols": 10,
                "nrows": 10,
                "xllcorner": 0.0,
                "yllcorner": 0.0,
                "cellsize": 100.0,
                "nodata_value": -9999.0,
            },
        )


def test_file_jobs_create_hydrology_and_watershed_outputs(tmp_path):
    source = tmp_path / "dem.tif"
    _dem(source)

    filled = fill_dem_file(source, tmp_path / "filled.tif")
    terrain = terrain_files(
        filled["filled_dem"], tmp_path / "slope.tif", tmp_path / "aspect.tif"
    )
    hydrology = hydrology_files(
        filled["filled_dem"],
        accumulation_path=tmp_path / "accumulation.tif",
        area_path=tmp_path / "area.tif",
        direction_path=tmp_path / "direction.tif",
        channel_path=tmp_path / "channel.shp",
        threshold_cells=1,
    )
    watershed = delineate_outlet_file(
        filled["filled_dem"],
        950,
        50,
        tmp_path / "watershed.tif",
        tmp_path / "watershed.shp",
    )

    outputs = (
        *terrain.values(),
        *hydrology.values(),
        watershed["raster_path"],
        watershed["vector_path"],
    )
    assert all(os.path.isfile(path) for path in outputs)
    assert watershed["cell_center"] == (950.0, 50.0)
    assert watershed["catchment_area_m2"] > 0
    domain_dem = materialize_domain_dem_file(
        filled["filled_dem"],
        watershed["vector_path"],
        tmp_path / "outlet" / "dem.asc",
    )
    assert os.path.isfile(domain_dem)
    domain_dataset = gdal.Open(domain_dem)
    assert domain_dataset.RasterXSize <= 10
    assert domain_dataset.RasterYSize <= 10
    assert domain_dataset.GetRasterBand(1).GetNoDataValue() == -9999

    domains = delineate_domains_file(
        filled["filled_dem"],
        (
            (
                "outlet",
                950,
                50,
                tmp_path / "domain.tif",
                tmp_path / "domain.shp",
                1,
                tmp_path / "outlet-2" / "dem.asc",
            ),
        ),
        dem_domain=(
            2,
            tmp_path / "dem_domain.tif",
            tmp_path / "dem_domain.shp",
            tmp_path / "dem_extent" / "dem.asc",
        ),
        merged_path=tmp_path / "merged.shp",
    )
    assert os.path.isfile(domains["merged_path"])
    assert domains["outlets"]["outlet"]["catchment_area_m2"] > 0
    # Delineation only records where the domain DEM will go; Morphology Setup
    # writes it later on the common L0 grid so every domain shares the extent.
    assert domains["outlets"]["outlet"]["dem_path"].endswith("outlet-2/dem.asc")
    assert not os.path.exists(domains["outlets"]["outlet"]["dem_path"])
    assert domains["dem_domain_path"].endswith("dem_extent/dem.asc")
    assert not os.path.exists(domains["dem_domain_path"])


def _degree_dem(path, cellsize, cols=24, rows=12):
    """Write a lat/lon DEM whose cell size is a repeating decimal."""
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path), cols, rows, 1, gdal.GDT_Float32
    )
    dataset.SetGeoTransform((79.0, cellsize, 0, 31.0, 0, -cellsize))
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(4326)
    dataset.SetProjection(crs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(
        np.arange(cols * rows, dtype="float32").reshape(rows, cols)
    )
    dataset.GetRasterBand(1).SetNoDataValue(-9999)
    dataset = None


def test_repeating_degree_cell_size_stays_aligned_but_a_floored_one_does_not(tmp_path):
    """Regression: flooring 1/1200 deg to 8 places drifts off the source grid."""
    from pymhm.grid_resolution import aligned_l0_l2_headers

    exact = 1.0 / 1200.0
    source = tmp_path / "dem.tif"
    _degree_dem(source, exact)
    anchor = {"xllcorner": 79.0, "yllcorner": 31.0 - 12 * exact, "unit": "deg"}
    bounds = (79.0, 79.0 + 24 * exact, 31.0 - 12 * exact, 31.0)

    l0_header, l2_header = aligned_l0_l2_headers(
        bounds, {**anchor, "cellsize": exact}, multiplier=12, unit="deg"
    )
    assert l0_header["ncols"] == l2_header["ncols"] * 12
    output = tmp_path / "crop.tif"
    crop_aligned_l0_raster(source, output, l0_header, reference_path=source)
    np.testing.assert_array_equal(
        gdal.Open(str(output)).ReadAsArray(),
        gdal.Open(str(source)).ReadAsArray(),
    )

    # The same grid built from an 8-decimal cell size no longer lands on the
    # source cells, and the window copy refuses it instead of resampling.
    floored_l0, _floored_l2 = aligned_l0_l2_headers(
        bounds, {**anchor, "cellsize": 0.00083333}, multiplier=12, unit="deg"
    )
    with pytest.raises(ValueError, match="aligned"):
        crop_aligned_l0_raster(
            source, tmp_path / "bad.tif", floored_l0, reference_path=source
        )

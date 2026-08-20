"""Domain DEMs must share the common L0 grid with every other model input."""
from __future__ import annotations

import os

import numpy as np
import pytest
from osgeo import gdal, ogr, osr

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pymhm import standalone_qgis  # noqa: E402

standalone_qgis.install(force=True)

from pymhm.Morphology.file_tasks import write_domain_dem_ascii  # noqa: E402
from pymhm.Morphology.latlon.ascii_pad import read_ascii_header  # noqa: E402


L0 = {
    "ncols": 6, "nrows": 5, "xllcorner": 0.0, "yllcorner": 0.0,
    "cellsize": 100.0, "nodata_value": -9999.0,
}


def _cropped_dem(path):
    """A DEM already cropped to the common L0 extent."""
    dataset = gdal.GetDriverByName("GTiff").Create(
        str(path), L0["ncols"], L0["nrows"], 1, gdal.GDT_Float32
    )
    dataset.SetGeoTransform((0.0, 100.0, 0, 500.0, 0, -100.0))
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(32632)
    dataset.SetProjection(crs.ExportToWkt())
    values = np.arange(
        L0["ncols"] * L0["nrows"], dtype="float32"
    ).reshape(L0["nrows"], L0["ncols"]) + 1
    dataset.GetRasterBand(1).WriteArray(values)
    dataset.GetRasterBand(1).SetNoDataValue(-9999)
    dataset = None
    return values


def _polygon(path, xmin, xmax, ymin, ymax):
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(32632)
    source = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(str(path))
    layer = source.CreateLayer("domain", crs, ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in ((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)):
        ring.AddPoint(x, y)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(polygon)
    layer.CreateFeature(feature)
    source = None


def _rows(path):
    lines = path.read_text(encoding="utf-8").splitlines()[6:]
    return [line.split() for line in lines if line.strip()]


def test_the_domain_dem_keeps_the_common_extent_not_the_domain_bounding_box(tmp_path):
    """This is the point of the change: same matrix size for every domain."""
    dem = tmp_path / "1_dem_filled_crop.tif"
    polygon = tmp_path / "domain.shp"
    output = tmp_path / "data" / "outlet_1" / "dem.asc"
    values = _cropped_dem(dem)
    # x 200-400 covers column centres 250/350 (columns 2-3); with ymax 500,
    # y 200-400 covers row centres 350/250 (rows 1-2).
    _polygon(polygon, 200.0, 400.0, 200.0, 400.0)

    write_domain_dem_ascii(dem, output, L0, polygon, reference_path=dem)

    header = read_ascii_header(output)
    assert (int(header["ncols"]), int(header["nrows"])) == (6, 5)
    assert header["xllcorner"] == 0.0 and header["yllcorner"] == 0.0
    assert header["cellsize"] == 100.0

    rows = _rows(output)
    assert len(rows) == 5 and all(len(row) == 6 for row in rows)
    kept = np.zeros((5, 6), dtype=bool)
    kept[1:3, 2:4] = True
    written = np.array([[float(cell) for cell in row] for row in rows])
    np.testing.assert_array_equal(written[kept], values[kept])
    assert np.all(written[~kept] == -9999)


def test_two_domains_produce_identical_matrix_sizes(tmp_path):
    dem = tmp_path / "crop.tif"
    _cropped_dem(dem)
    headers = []
    for index, bounds in enumerate(
            ((0.0, 200.0, 300.0, 500.0), (400.0, 600.0, 0.0, 200.0))):
        polygon = tmp_path / f"d{index}.shp"
        _polygon(polygon, *bounds)
        output = tmp_path / f"d{index}" / "dem.asc"
        write_domain_dem_ascii(dem, output, L0, polygon, reference_path=dem)
        headers.append(read_ascii_header(output))
    # Different catchments, same grid -- which is what mHM needs.
    assert headers[0] == headers[1]


def test_a_domain_outside_the_extent_yields_an_all_nodata_grid(tmp_path):
    dem = tmp_path / "crop.tif"
    polygon = tmp_path / "far.shp"
    output = tmp_path / "dem.asc"
    _cropped_dem(dem)
    _polygon(polygon, 10000.0, 10200.0, 10000.0, 10200.0)

    write_domain_dem_ascii(dem, output, L0, polygon, reference_path=dem)

    rows = _rows(output)
    assert len(rows) == 5
    assert all(cell == "-9999" for row in rows for cell in row)


def test_a_failed_write_leaves_no_partial_file(tmp_path):
    dem = tmp_path / "crop.tif"
    output = tmp_path / "dem.asc"
    _cropped_dem(dem)

    with pytest.raises(RuntimeError, match="Could not open the mask layer"):
        write_domain_dem_ascii(dem, output, L0, tmp_path / "missing.shp")
    assert not output.exists()
    assert not list(tmp_path.glob(".*.tmp"))


def test_source_nodata_becomes_the_ascii_nodata(tmp_path):
    dem = tmp_path / "crop.tif"
    polygon = tmp_path / "all.shp"
    output = tmp_path / "dem.asc"
    _cropped_dem(dem)
    # Punch a nodata hole in the source.
    dataset = gdal.Open(str(dem), gdal.GA_Update)
    band = dataset.GetRasterBand(1)
    values = band.ReadAsArray()
    values[0, 0] = -9999
    band.WriteArray(values)
    dataset = None
    _polygon(polygon, 0.0, 600.0, 0.0, 500.0)

    write_domain_dem_ascii(dem, output, L0, polygon, reference_path=dem)

    assert _rows(output)[0][0] == "-9999"


def _domain_state(project, outlets, dem_domain=False):
    from pymhm.Morphology.watershed.domain_state import save_state

    save_state(project, {
        "definition_mode": "snapped_pour_points",
        "pour_points_source": "outlets.shp",
        "outlet_id_field": "id",
        "outlet_order": list(outlets),
        "dem_domain": dem_domain,
        "outlets": {
            outlet: {"is_domain": True, "vector_path": f"vec/{outlet}.shp"}
            for outlet in outlets
        },
    })


def test_the_plan_lists_every_active_domain_with_its_polygon(tmp_path):
    from pymhm.Morphology.layers.domain_dem import domain_dem_plan

    _domain_state(tmp_path, ("001", "002"), dem_domain=True)
    plan = domain_dem_plan(tmp_path)

    assert [entry["name"] for entry in plan] == ["001", "002", "dem_extent"]
    assert [entry["domain_id"] for entry in plan] == [1, 2, 3]
    # Each entry names where its DEM will be written, on the shared grid.
    for entry in plan:
        assert entry["dem_path"].endswith(f"data/{entry['name']}/dem.asc")
        assert entry["polygon"]
    # The valid-DEM domain uses the delineator polygon, not an outlet vector.
    assert plan[-1]["polygon"].endswith("4_watershed_DEM.shp")


def test_writing_domain_dems_refuses_a_missing_polygon(tmp_path):
    from pymhm.Morphology.layers.domain_dem import write_domain_dems

    dem = tmp_path / "crop.tif"
    _cropped_dem(dem)
    _domain_state(tmp_path, ("001",))

    with pytest.raises(FileNotFoundError, match="no delineated polygon"):
        write_domain_dems(tmp_path, str(dem), L0)


def test_writing_domain_dems_refuses_a_missing_cropped_dem(tmp_path):
    from pymhm.Morphology.layers.domain_dem import write_domain_dems

    _domain_state(tmp_path, ("001",))
    with pytest.raises(FileNotFoundError, match="cropped L0 DEM is required"):
        write_domain_dems(tmp_path, str(tmp_path / "absent.tif"), L0)


def test_no_active_domains_is_a_no_op(tmp_path):
    from pymhm.Morphology.layers.domain_dem import write_domain_dems

    assert write_domain_dems(tmp_path, str(tmp_path / "absent.tif"), L0) == []


def test_the_state_records_the_domain_plan(tmp_path):
    from pymhm.Morphology.core.processing_state import ProcessingStateMixin
    from pymhm.state_cache import load_state

    class _State(ProcessingStateMixin):
        def __init__(self, project):
            self.dialog = type("D", (), {"project_folder": str(project)})()
            self.processing_state_filename = "pymhm_processing_state.json"
            self.processing_state = {"version": 1, "outputs": {}, "workflows": {}}
            self.log_message = lambda _m: None

    state = _State(tmp_path)
    plan = [{
        "domain_id": 1, "outlet_id": "001", "name": "001",
        "polygon": "/p/vec/001.shp", "directory": "/p/data/001",
        "dem_path": "/p/data/001/dem.asc",
    }]
    state.save_domain_plan(plan)

    assert state.saved_domain_plan() == plan
    assert load_state(tmp_path)["domains"] == plan

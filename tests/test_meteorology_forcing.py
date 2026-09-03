from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from mhm_qgis.core.meteorology.forcing import (
    ERA5LAND,
    MHM_READY,
    MeteoFolderSpec,
    TargetGrid,
    inspect_meteo_folder,
    process_meteo_inputs,
    resolution_in_crs,
)


LAT = np.array([51.0, 50.0])
LON = np.array([10.0, 11.0])


def _write_nc(
        path: Path,
        variable: str,
        times=("2001-01-01T12:00:00",),
        lat=LAT,
        lon=LON,
        values=None,
        units=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    times = np.asarray(times, dtype="datetime64[ns]")
    if values is None:
        values = np.ones((len(times), len(lat), len(lon)), dtype="float64")
    ds = xr.Dataset(
        {variable: (("time", "latitude", "longitude"), values)},
        coords={"time": times, "latitude": lat, "longitude": lon},
    )
    if units:
        ds[variable].attrs["units"] = units
    ds.to_netcdf(path, engine="scipy")


def _write_projected_nc(path: Path, variable: str, value: float):
    from pyproj import CRS

    path.parent.mkdir(parents=True, exist_ok=True)
    ds = xr.Dataset(
        {
            variable: (
                ("time", "yc", "xc"),
                np.full((1, 2, 2), value, dtype="float64"),
                {"grid_mapping": "crs"},
            ),
            "crs": ((), 0, {"spatial_ref": CRS.from_epsg(3857).to_wkt()}),
        },
        coords={
            "time": np.asarray(["2001-01-01"], dtype="datetime64[ns]"),
            "yc": [1500.0, 500.0],
            "xc": [500.0, 1500.0],
        },
    )
    ds.to_netcdf(path, engine="scipy")


def test_mhm_ready_folder_contracts_and_resolution(tmp_path):
    pre = tmp_path / "pre"
    temp = tmp_path / "temp"
    pet = tmp_path / "pet"
    _write_nc(pre / "rain.nc", "pre")
    for variable in ("tavg", "tmin", "tmax"):
        _write_nc(temp / f"{variable}.nc", variable)
    _write_nc(pet / "pet.nc", "pet")

    metadata = inspect_meteo_folder(pre, "precipitation", "mHM_ready")
    assert metadata.source == MHM_READY
    assert metadata.shape == (2, 2)
    assert metadata.resolution == pytest.approx(1.0)
    assert metadata.unit == "deg"
    assert inspect_meteo_folder(temp, "temperature", MHM_READY).shape == (2, 2)
    assert inspect_meteo_folder(pet, "pet", MHM_READY).files == (
        pet / "pet.nc",)

    projected = resolution_in_crs(metadata, "EPSG:3857")
    assert projected.unit == "m"
    assert 100_000 < projected.resolution < 200_000

    _write_nc(pre / "extra.nc", "pre")
    with pytest.raises(ValueError, match="exactly one"):
        inspect_meteo_folder(pre, "precipitation", MHM_READY)


def test_mhm_ready_temperature_requires_exact_files(tmp_path):
    folder = tmp_path / "temperature"
    _write_nc(folder / "tavg.nc", "tavg")
    _write_nc(folder / "tmin.nc", "tmin")
    _write_nc(folder / "unexpected.nc", "temperature")

    with pytest.raises(ValueError, match="folder must contain exactly") as error:
        inspect_meteo_folder(folder, "temperature", MHM_READY)
    assert "unexpected.nc" not in str(error.value)


def test_era5land_rejects_incompatible_grids(tmp_path):
    folder = tmp_path / "pre"
    _write_nc(folder / "a.nc", "tp")
    _write_nc(folder / "b.nc", "tp", times=("2001-01-02T12:00:00",),
              lon=np.array([10.0, 12.0]))

    with pytest.raises(ValueError, match="identical spatial coordinates"):
        inspect_meteo_folder(folder, "precipitation", ERA5LAND)


def test_processes_separate_era5land_folders_and_optional_pet(
    tmp_path,
    monkeypatch,
):
    pre = tmp_path / "raw" / "pre"
    temp = tmp_path / "raw" / "temp"
    pet = tmp_path / "raw" / "pet"
    for index, day in enumerate((1, 2), start=1):
        times = tuple(
            f"2001-01-{day:02d}T{hour:02d}:00:00"
            for hour in (6, 12, 18))
        pre_values = np.asarray([
            [[0.001, -2e-6], [-0.0, 1e-12]],
            [[0.002, -1e-6], [-0.0, 1e-12]],
            [[0.0015, -1.5e-6], [-0.0, 1e-12]],
        ])
        temp_values = np.stack([
            np.full((2, 2), value)
            for value in (273.15, 275.15, 277.15)])
        pet_values = np.stack([
            np.full((2, 2), value)
            for value in (-0.001, -0.002, -0.003)])
        _write_nc(
            pre / f"{index}.nc", "tp", times=times,
            values=pre_values, units="m")
        _write_nc(
            temp / f"{index}.nc", "t2m", times=times,
            values=temp_values, units="K")
        _write_nc(
            pet / f"{index}.nc", "pev", times=times,
            values=pet_values, units="m")

    header = {
        "ncols": 2,
        "nrows": 2,
        "xllcorner": 9.5,
        "yllcorner": 49.5,
        "cellsize": 1.0,
        "nodata_value": -9999.0,
    }
    result = process_meteo_inputs(
        MeteoFolderSpec("precipitation", pre, ERA5LAND),
        MeteoFolderSpec("temperature", temp, ERA5LAND),
        tmp_path / "out",
        TargetGrid(lon=LON, lat=LAT, header=header),
        pet=MeteoFolderSpec("pet", pet, ERA5LAND),
    )

    assert set(result.outputs) == {"pre", "tavg", "tmin", "tmax", "pet"}
    for variable, output in result.outputs.items():
        assert output == tmp_path / "out" / variable / f"{variable}.nc"
        assert output.exists()
        assert result.headers[variable].read_text().startswith("ncols")
        with xr.open_dataset(output) as ds:
            assert ds[variable].dims == ("time", "lat", "lon")
            assert ds[variable].shape == (2, 2, 2)

    with xr.open_dataset(result.outputs["pre"]) as ds:
        assert float(ds["pre"].isel(time=0, lat=0, lon=0)) == pytest.approx(2.0)
        assert float(ds["pre"].isel(time=0, lat=0, lon=1)) == 0.0
        assert float(ds["pre"].isel(time=0, lat=1, lon=0)) == 0.0
        assert not np.signbit(ds["pre"].values[ds["pre"].values == 0]).any()
        assert float(ds["pre"].isel(time=0, lat=1, lon=1)) == pytest.approx(
            1e-9)
        assert ds["pre"].attrs["valid_min"] == 0.0
        assert "below 0 mm clipped to 0" in ds.attrs["history"]
    with xr.open_dataset(result.outputs["tavg"]) as ds:
        assert float(ds["tavg"].isel(time=0, lat=0, lon=0)) == pytest.approx(2.0)
    with xr.open_dataset(result.outputs["tmin"]) as ds:
        assert float(ds["tmin"].isel(time=0, lat=0, lon=0)) == pytest.approx(0.0)
    with xr.open_dataset(result.outputs["tmax"]) as ds:
        assert float(ds["tmax"].isel(time=0, lat=0, lon=0)) == pytest.approx(4.0)
    with xr.open_dataset(result.outputs["pet"]) as ds:
        assert float(ds["pet"].isel(time=0, lat=0, lon=0)) == pytest.approx(3.0)

    without_pet = process_meteo_inputs(
        MeteoFolderSpec("precipitation", pre, ERA5LAND),
        MeteoFolderSpec("temperature", temp, ERA5LAND),
        tmp_path / "out",
        TargetGrid(lon=LON, lat=LAT, header=header),
    )
    assert set(without_pet.outputs) == {"pre", "tavg", "tmin", "tmax"}
    assert not (tmp_path / "out" / "pet" / "pet.nc").exists()
    assert not (tmp_path / "out" / "pet" / "header.txt").exists()

    from mhm_qgis.core.meteorology import forcing

    previous_pre = without_pet.outputs["pre"].read_bytes()
    original_writer = forcing.write_netcdf

    def fail_after_staging_some_outputs(ds, variable, output):
        if variable == "tmin":
            raise RuntimeError("synthetic staging failure")
        return original_writer(ds, variable, output)

    monkeypatch.setattr(
        forcing,
        "write_netcdf",
        fail_after_staging_some_outputs,
    )
    with pytest.raises(RuntimeError, match="synthetic staging failure"):
        process_meteo_inputs(
            MeteoFolderSpec("precipitation", pre, ERA5LAND),
            MeteoFolderSpec("temperature", temp, ERA5LAND),
            tmp_path / "out",
            TargetGrid(lon=LON, lat=LAT, header=header),
        )
    assert without_pet.outputs["pre"].read_bytes() == previous_pre
    assert not list((tmp_path / "out").glob(".mhm_qgis-meteo-*"))


def test_mhm_ready_only_clips_precipitation(tmp_path):
    pre = tmp_path / "ready" / "pre"
    temp = tmp_path / "ready" / "temp"
    pet = tmp_path / "ready" / "pet"
    _write_nc(
        pre / "pre.nc",
        "pre",
        values=np.asarray([[[-1e-9, -0.0], [np.nan, 1e-12]]]),
        units="mm",
    )
    for variable, value in (("tavg", -5.0), ("tmin", -10.0), ("tmax", -1.0)):
        _write_nc(
            temp / f"{variable}.nc",
            variable,
            values=np.full((1, 2, 2), value),
            units="Celsius",
        )
    _write_nc(
        pet / "pet.nc",
        "pet",
        values=np.full((1, 2, 2), -2.0),
        units="mm",
    )
    header = {
        "ncols": 2,
        "nrows": 2,
        "xllcorner": 9.5,
        "yllcorner": 49.5,
        "cellsize": 1.0,
        "nodata_value": -9999.0,
    }

    result = process_meteo_inputs(
        MeteoFolderSpec("precipitation", pre, MHM_READY),
        MeteoFolderSpec("temperature", temp, MHM_READY),
        tmp_path / "out",
        TargetGrid(lon=LON, lat=LAT, header=header),
        pet=MeteoFolderSpec("pet", pet, MHM_READY),
    )

    with xr.open_dataset(result.outputs["pre"]) as ds:
        values = ds["pre"].values
        assert values[0, 0, 0] == 0.0
        assert values[0, 0, 1] == 0.0
        assert not np.signbit(values[values == 0]).any()
        assert np.isnan(values[0, 1, 0])
        assert values[0, 1, 1] == pytest.approx(1e-12)
        assert ds["pre"].attrs["valid_min"] == 0.0
    with xr.open_dataset(result.outputs["tmin"]) as ds:
        assert np.all(ds["tmin"].values == -10.0)
        assert "valid_min" not in ds["tmin"].attrs
    with xr.open_dataset(result.outputs["pet"]) as ds:
        assert np.all(ds["pet"].values == -2.0)
        assert "valid_min" not in ds["pet"].attrs


def test_processes_projected_mhm_ready_files_on_exact_header_grid(tmp_path):
    from pyproj import Transformer

    pre = tmp_path / "ready" / "pre"
    temp = tmp_path / "ready" / "temp"
    _write_projected_nc(pre / "pre.nc", "pre", 4.0)
    for index, variable in enumerate(("tavg", "tmin", "tmax")):
        _write_projected_nc(temp / f"{variable}.nc", variable, float(index))

    header = {
        "ncols": 2,
        "nrows": 2,
        "xllcorner": 0.0,
        "yllcorner": 0.0,
        "cellsize": 1000.0,
        "nodata_value": -9999.0,
    }
    transform = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    x_mesh, y_mesh = np.meshgrid(
        np.asarray([500.0, 1500.0]),
        np.asarray([1500.0, 500.0]),
    )
    sample_lon, sample_lat = transform.transform(x_mesh, y_mesh)
    lon = sample_lon[0, :]
    lat = sample_lat[:, 0]
    result = process_meteo_inputs(
        MeteoFolderSpec("precipitation", pre, MHM_READY),
        MeteoFolderSpec("temperature", temp, MHM_READY),
        tmp_path / "out",
        TargetGrid(
            lon=lon,
            lat=lat,
            header=header,
            crs="EPSG:3857",
            sample_lon=sample_lon,
            sample_lat=sample_lat,
        ),
    )

    with xr.open_dataset(result.outputs["pre"]) as ds:
        assert ds["pre"].shape == (1, 2, 2)
        assert np.all(ds["pre"].values == 4.0)
        assert ds["lon2d"].shape == (2, 2)
        assert ds["lat2d"].shape == (2, 2)
        assert "crs" in ds
        assert ds["pre"].attrs["grid_mapping"] == "crs"

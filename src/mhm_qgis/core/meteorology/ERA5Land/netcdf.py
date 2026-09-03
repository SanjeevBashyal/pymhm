from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, List, Union

import numpy as np
import xarray as xr
import geopandas as gpd
import rioxarray  # noqa: F401 - needed for rio accessor


class NetCDFProcessor:
    """Utilities for working with NetCDF datasets and vector AOIs."""

    def __init__(self, nc_path: Optional[Path] = None, ds: Optional[xr.Dataset] = None) -> None:
        if (nc_path is None) == (ds is None):
            raise ValueError("Provide exactly one of nc_path or ds")
        self.nc_path = Path(nc_path) if nc_path is not None else None
        self.ds = xr.open_dataset(
            self.nc_path, engine="netcdf4") if ds is None else ds
        # Ensure CRS for rioxarray
        self.ds = self.ds.rio.write_crs("epsg:4326", inplace=True)

    def detect_lat_lon_names(self) -> Tuple[str, str]:
        lat_candidates = ["lat", "latitude", "y"]
        lon_candidates = ["lon", "longitude", "x"]
        lat_name = next(
            (n for n in lat_candidates if n in self.ds.coords), None)
        lon_name = next(
            (n for n in lon_candidates if n in self.ds.coords), None)
        if lat_name is None or lon_name is None:
            raise KeyError("Latitude/longitude coordinates not found.")
        return lat_name, lon_name

    def subset_to_vector_bounds(self, vector_path: Path) -> xr.Dataset:
        gdf = gpd.read_file(vector_path)
        if gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)
        minx, miny, maxx, maxy = gdf.total_bounds
        lat_name, lon_name = self.detect_lat_lon_names()
        lon_vals = self.ds[lon_name].values
        if lon_vals.min() >= 0 and lon_vals.max() <= 360 and minx < 0:
            minx += 360
            maxx += 360
        lat_vals = self.ds[lat_name].values
        lat_desc = lat_vals[0] > lat_vals[-1]
        lat_slice = slice(maxy, miny) if lat_desc else slice(miny, maxy)
        self.ds = self.ds.sel(
            {lon_name: slice(minx, maxx), lat_name: lat_slice})
        return self.ds

    def clip_to_vector(self, vector_path: Path) -> xr.Dataset:
        gdf = gpd.read_file(vector_path)
        if gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)
        try:
            lat_name, lon_name = self.detect_lat_lon_names()
            try:
                self.ds = self.ds.rio.set_spatial_dims(
                    x_dim=lon_name, y_dim=lat_name, inplace=False)
            except TypeError:
                self.ds = self.ds.rio.set_spatial_dims(
                    lon_name, lat_name, inplace=False)
        except Exception:
            pass
        return self.ds.rio.clip(gdf.geometry, gdf.crs, drop=True, all_touched=True)

    def extract_point(self, lat: float, lon: float, var: Optional[str] = None):
        lat_name, lon_name = self.detect_lat_lon_names()
        lon_vals = self.ds[lon_name].values
        if lon_vals.min() >= 0 and lon_vals.max() <= 360 and lon < 0:
            lon = lon + 360
        selected = self.ds.sel(
            {lat_name: lat, lon_name: lon}, method="nearest")
        if var:
            if var not in selected:
                raise KeyError(f"Variable '{var}' not found")
            return selected[var]
        return selected

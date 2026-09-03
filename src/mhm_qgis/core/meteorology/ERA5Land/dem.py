from __future__ import annotations

from typing import Optional, Iterable, List

import geopandas as gpd
import rasterio
from shapely.geometry import Point


class DEMSampler:
    def __init__(self, dem_path: str) -> None:
        self.dem_path = dem_path

    def sample(self, lat: float, lon: float) -> Optional[float]:
        try:
            point_gdf = gpd.GeoDataFrame(geometry=[Point(lon, lat)], crs="EPSG:4326")
            with rasterio.open(self.dem_path) as src:
                point_proj = point_gdf.to_crs(src.crs)
                coords = [(p.x, p.y) for p in point_proj.geometry]
                elevation_values = list(src.sample(coords))
                return float(elevation_values[0][0])
        except Exception:
            return None

    def sample_many(self, lats: Iterable[float], lons: Iterable[float]) -> List[Optional[float]]:
        """Sample elevations for many points efficiently using one raster open.

        Returns a list aligned with the input iterables.
        """
        try:
            points = [Point(lon, lat) for lat, lon in zip(lats, lons)]
            gdf = gpd.GeoDataFrame(geometry=points, crs="EPSG:4326")
            with rasterio.open(self.dem_path) as src:
                gdf_proj = gdf.to_crs(src.crs)
                coords = [(p.x, p.y) for p in gdf_proj.geometry]
                values = list(src.sample(coords))
                return [float(v[0]) if v is not None else None for v in values]
        except Exception:
            # On any failure, fall back to per-point sampling (slower but safer)
            return [self.sample(lat, lon) for lat, lon in zip(lats, lons)]



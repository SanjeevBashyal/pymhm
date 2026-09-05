"""CRS units for local distance calculations, independent of QGIS."""


def local_distance_crs(crs, x, y):
    """Return a distance CRS and metres per unit at a project coordinate.

    Geographic searches use an outlet-centred equidistant projection: a metre
    buffer transforms back to degrees at this latitude, not a fixed degree radius.
    """
    from pyproj import CRS, Transformer

    source = CRS.from_user_input(crs)
    if source.is_geographic:
        lon, lat = Transformer.from_crs(source, "EPSG:4326", always_xy=True).transform(x, y)
        local = CRS.from_proj4(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m")
        return local.to_wkt(), 1.0
    if not source.is_projected:
        raise ValueError("Snapping requires a projected or geographic project CRS.")
    return source.to_wkt(), source.axis_info[0].unit_conversion_factor

"""Variable aliases used when reading local ERA5-Land NetCDF files."""

VARIABLE_NAME_MAPPINGS = {
    "precipitation": [
        "tp",
        "Total_precipitation_surface_1_Hour_Accumulation",
        "total_precipitation",
    ],
    "temperature": [
        "t2m",
        "Temperature_height_above_ground",
        "2m_temperature",
    ],
    "potential_evaporation": ["pev", "potential_evaporation"],
    "dewpoint": [
        "d2m",
        "Dewpoint_temperature_height_above_ground",
        "2m_dewpoint_temperature",
    ],
    "solar_radiation": [
        "ssrd",
        "surface_solar_radiation_downwards",
        "Downward_Short-Wave_Radiation_Flux_surface",
    ],
    "u_wind": [
        "u10",
        "10m_u_component_of_wind",
        "u-component_of_wind_height_above_ground",
    ],
    "v_wind": [
        "v10",
        "10m_v_component_of_wind",
        "v-component_of_wind_height_above_ground",
    ],
}

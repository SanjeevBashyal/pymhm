"""Tests for the project-local namelist handoff."""

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

from pymhm.Morphology.watershed.domain_state import save_state
from pymhm.nml_settings import load_settings, sync_domain_settings, update_section
from pymhm.vSpecific import build_initial_values


class Dialog:
    def __init__(self, project_folder):
        self.project_folder = str(project_folder)

    def current_l1_resolution(self):
        return 1000

    def current_l11_resolution(self):
        return 2000

    def current_grid_unit(self):
        return "m"


def test_sections_are_merged_without_overwriting_future_inputs(tmp_path):
    update_section(tmp_path, "land_cover", {"variable": "land_cover"})
    update_section(tmp_path, "soil", {"variable": "soil_class"})

    state = load_settings(tmp_path)
    assert state["land_cover"]["variable"] == "land_cover"
    assert state["soil"]["variable"] == "soil_class"


def test_domain_gauges_project_to_v5_namelist_values(tmp_path):
    save_state(
        tmp_path,
        {
            "definition_mode": "snapped_pour_points",
            "pour_points_source": "pour_points.shp",
            "outlet_id_field": "station_id",
            "outlet_order": ["001"],
            "dem_domain": False,
            "outlets": {
                "001": {
                    "is_domain": True,
                    "is_gauged": True,
                    "gauge_id": 1,
                    "gauge_filename": "001.txt",
                    "gauge_path": "data/observation/streamflow/001.txt",
                    "domain_ids": [1],
                }
            },
        },
    )
    sync_domain_settings(tmp_path)

    state = load_settings(tmp_path)
    assert state["domains"] == [
        {"domain_id": 1, "is_dem_domain": False, "outlet_id": "001"}
    ]
    assert state["gauges"][0]["gauge_filename"] == "001.txt"

    gauges = build_initial_values("5.13", Dialog(tmp_path))["main"][
        "evaluation_gauges"
    ]
    assert gauges["nGaugesTotal"] == 1
    assert gauges["NoGauges_domain"] == [1]
    assert gauges["Gauge_id"][0][0] == 1
    assert gauges["gauge_filename"][0][0] == "001.txt"

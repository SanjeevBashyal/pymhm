"""Tests for the compact morphology display registry."""

from pathlib import Path

from pymhm.morphology_display import land_cover_periods, resolve_display_output
from pymhm.nml_settings import update_section
from pymhm.project_layout import ensure_project_structure, morph_folder


def test_historical_land_cover_resolves_the_selected_period(tmp_path):
    ensure_project_structure(tmp_path, "5.13")
    morph = Path(morph_folder(tmp_path))
    first = morph / "lc_1990_1999.asc"
    second = morph / "lc_2000_2010.asc"
    first.write_text("first", encoding="utf-8")
    second.write_text("second", encoding="utf-8")
    update_section(
        tmp_path,
        "land_cover",
        {
            "scenes": [
                {
                    "start_year": 1990,
                    "end_year": 1999,
                    "output_path": "data/static/morph/lc_1990_1999.asc",
                },
                {
                    "start_year": 2000,
                    "end_year": 2010,
                    "output_path": "data/static/morph/lc_2000_2010.asc",
                },
            ]
        },
    )

    assert len(land_cover_periods(tmp_path)) == 2
    assert resolve_display_output(
        tmp_path, "land_cover", year=2005
    ).path == second

"""Tests for the plugin-owned project workspace layout."""
from __future__ import annotations

from pathlib import Path

from pymhm.project_layout import (
    data_folder,
    data_raw_folder,
    ensure_project_structure,
    geometry_folder,
    lai_folder,
    meteo_folder,
    morph_folder,
    output_folder,
    raw_meteo_folder,
    relative_project_path,
    restart_folder,
    static_folder,
    streamflow_observation_folder,
    workspace_folder,
    z_temp_folder,
)


def test_all_generated_paths_are_inside_plugin_workspace(tmp_path: Path) -> None:
    project = tmp_path / "project"
    workspace = project / "mhm-plugin"

    assert Path(workspace_folder(project)) == workspace
    assert Path(z_temp_folder(project)) == workspace / "Z Temp"
    assert Path(geometry_folder(project)) == workspace / "Z Temp" / "Geometry"
    assert Path(data_folder(project)) == workspace / "data"
    assert Path(data_raw_folder(project)) == workspace / "data_raw"
    assert Path(static_folder(project)) == workspace / "data" / "static"
    assert Path(morph_folder(project)) == workspace / "data" / "static" / "morph"
    assert Path(meteo_folder(project)) == workspace / "data" / "meteo"
    assert Path(raw_meteo_folder(project)) == workspace / "data_raw" / "meteo"
    assert Path(lai_folder(project)) == workspace / "data" / "lai"
    assert (
        Path(streamflow_observation_folder(project))
        == workspace / "data" / "observation" / "streamflow"
    )
    assert Path(output_folder(project)) == workspace / "output"
    assert Path(restart_folder(project)) == workspace / "restart"


def test_ensure_project_structure_does_not_create_outputs_at_outer_root(
    tmp_path: Path,
) -> None:
    project = tmp_path / "project"

    ensure_project_structure(project, "6.0")

    workspace = project / "mhm-plugin"
    assert (workspace / "Z Temp" / "Geometry").is_dir()
    assert (workspace / "data" / "static" / "morph").is_dir()
    assert (workspace / "data" / "meteo").is_dir()
    assert (workspace / "data" / "lai").is_dir()
    assert (workspace / "output").is_dir()
    assert (workspace / "restart").is_dir()
    for name in ("Z Temp", "data", "data_raw", "output", "restart"):
        assert not (project / name).exists()


def test_generated_path_is_relative_to_outer_project(tmp_path: Path) -> None:
    project = tmp_path / "project"
    output = Path(morph_folder(project)) / "dem.asc"

    assert (
        relative_project_path(project, output)
        == "mhm-plugin/data/static/morph/dem.asc"
    )

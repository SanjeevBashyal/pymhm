"""Focused tests for morphology workflow boundaries."""

from pathlib import Path

from pymhm import standalone_qgis

standalone_qgis.install(force=True)

# isort: off
from pymhm.Morphology.orchestration.execute_all import ExecuteAllMixin  # noqa: E402
from pymhm.project_layout import geometry_folder, morph_folder  # noqa: E402
# isort: on


class _WorkflowHarness:
    def __init__(self, project):
        self.dialog = type("Dialog", (), {"project_folder": str(project)})()
        self.calls = []
        self.statuses = []
        self.skip_loading = False
        self.categorical_ready_outputs = {}
        self.filled_dem_path = ""
        self.flow_accumulation_path = ""
        self.flow_direction_path = ""
        self.channel_network_vector_path = ""
        self.snapped_points_path = ""
        self.merged_watershed_path = ""
        self.latlon_result = True

    def _touch(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return str(path)

    def log_message(self, _message):
        pass

    def mark_workflow_status(self, workflow, status, message=""):
        self.statuses.append((workflow, status, message))

    def check_prerequisites(self, needs_pour_points=False):
        return True

    def without_layer_loading(self, callback):
        callback()

    def fill_dem(self):
        self.calls.append("fill")
        self.filled_dem_path = self._touch(
            Path(geometry_folder(self.dialog.project_folder))
            / "1_dem_filled.tif"
        )

    def process_slope(self):
        self.calls.append("slope")

    def process_aspect(self):
        self.calls.append("aspect")

    def process_land_use(self):
        self.calls.append("land_cover")
        return True

    def process_soil(self, write_classdefinition=True):
        self.calls.append("soil")
        geometry = Path(geometry_folder(self.dialog.project_folder))
        morph = Path(morph_folder(self.dialog.project_folder))
        self._touch(geometry / "3_soil.tif")
        self._touch(morph / "soil_classdefinition.txt")
        return True

    def process_geology(self, write_classdefinition=True):
        self.calls.append("geology")
        geometry = Path(geometry_folder(self.dialog.project_folder))
        morph = Path(morph_folder(self.dialog.project_folder))
        self._touch(geometry / "3_geology_processed.tif")
        self._touch(geometry / "geology_class_metadata.json")
        self._touch(morph / "geology_classdefinition.txt")
        return True

    def _categorical_mode(self, _kind):
        return "lookup_table"

    def process_flow_accumulation(self):
        self.calls.append("flow_accumulation")
        self.flow_accumulation_path = self._touch(
            Path(geometry_folder(self.dialog.project_folder))
            / "2_flow_accumulation.tif"
        )

    def process_flow_direction(self):
        self.calls.append("flow_direction")
        self.flow_direction_path = self._touch(
            Path(geometry_folder(self.dialog.project_folder))
            / "2_flow_direction.tif"
        )

    def process_channel_network(self):
        self.calls.append("channel_network")
        self.channel_network_vector_path = self._touch(
            Path(geometry_folder(self.dialog.project_folder))
            / "3_channel_network.gpkg"
        )

    def snap_points(self):
        self.calls.append("snap")
        self.snapped_points_path = self._touch(
            Path(geometry_folder(self.dialog.project_folder))
            / "3_snapped_points.gpkg"
        )

    def process_gauge_position(self):
        self.calls.append("gauge")

    def _restore_existing_path(self, _attribute, *_filenames):
        return None

    def delineate_watershed(self):
        self.calls.append("watershed")

    def crop_all_layers(self, show_error_dialog=True):
        self.calls.append("crop")
        return True

    def mask_all_layers(self, show_error_dialog=True):
        self.calls.append("mask")
        return True

    def process_lat_lon(self):
        self.calls.append("latlon")
        return self.latlon_result

    def write_all_layers(self, show_error_dialog=True):
        self.calls.append("write")
        return True


def test_execute_all_stops_before_finalization(tmp_path):
    workflow = _WorkflowHarness(tmp_path)

    assert ExecuteAllMixin.execute_all_processing(
        workflow,
        show_error_dialog=False,
    )
    assert workflow.calls == [
        "fill",
        "slope",
        "aspect",
        "land_cover",
        "soil",
        "geology",
        "flow_accumulation",
        "flow_direction",
        "channel_network",
        "snap",
        "gauge",
    ]
    assert workflow.statuses[-1][:2] == ("execute_all", "completed")


def test_combined_setup_runs_finalization_in_order(tmp_path):
    workflow = _WorkflowHarness(tmp_path)

    assert ExecuteAllMixin.execute_morph_setup_processing(
        workflow,
        show_error_dialog=False,
        workflow_key="meteo_morph_setup",
    )
    assert workflow.calls == ["crop", "mask", "latlon", "write"]
    assert workflow.statuses[-1][:2] == (
        "meteo_morph_setup",
        "completed",
    )


def test_combined_setup_stops_when_latlon_fails(tmp_path):
    workflow = _WorkflowHarness(tmp_path)
    workflow.latlon_result = False

    assert not ExecuteAllMixin.execute_morph_setup_processing(
        workflow,
        show_error_dialog=False,
        workflow_key="meteo_morph_setup",
    )
    assert workflow.calls == ["crop", "mask", "latlon"]
    assert workflow.statuses[-1][:2] == ("meteo_morph_setup", "failed")

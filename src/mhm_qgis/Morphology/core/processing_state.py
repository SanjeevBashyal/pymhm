# -*- coding: utf-8 -*-
"""Project-local registry for prepared morphology outputs."""
from ..common import (
    os,
    json,
    processing,
)
from ...project_layout import workspace_folder
from ...time_utils import utc_timestamp


class ProcessingStateMixin:
    """Project-local registry for prepared morphology outputs."""

    def processing_state_path(self):
        """Return the project-local processing state file path."""
        if not self.dialog.project_folder:
            return None
        return os.path.join(
            workspace_folder(self.dialog.project_folder),
            self.processing_state_filename
        )

    def load_processing_state(self):
        """Load the processing output registry for this project."""
        state_path = self.processing_state_path()
        if not state_path or not os.path.exists(state_path):
            self.processing_state = {
                "version": 1, "outputs": {}, "workflows": {}, "grid": {},
            }
            return

        try:
            with open(state_path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if not isinstance(state, dict):
                raise ValueError("Processing state is not a JSON object.")
            state.setdefault("version", 1)
            state.setdefault("outputs", {})
            state.setdefault("workflows", {})
            state.setdefault("grid", {})
            self.processing_state = state
        except Exception as e:
            self.log_message(f"WARNING: Could not read processing state: {e}")
            self.processing_state = {
                "version": 1, "outputs": {}, "workflows": {}, "grid": {},
            }

    OWNED_SECTIONS = ("version", "outputs", "workflows", "grid", "domains")

    def save_processing_state(self):
        """Write the sections this registry owns, preserving the others.

        `state_cache` and the meteorology state write fingerprints and cached
        inspections into the same file. Dumping the in-memory copy wholesale
        would erase whatever they added since it was loaded, silently undoing
        stage reuse, so only the owned sections are overlaid.
        """
        state_path = self.processing_state_path()
        if not state_path:
            return

        try:
            os.makedirs(os.path.dirname(state_path), exist_ok=True)
            merged = {}
            if os.path.exists(state_path):
                try:
                    with open(state_path, "r", encoding="utf-8") as state_file:
                        stored = json.load(state_file)
                    if isinstance(stored, dict):
                        merged = stored
                except Exception:
                    merged = {}
            for section in self.OWNED_SECTIONS:
                if section in self.processing_state:
                    merged[section] = self.processing_state[section]
            with open(state_path, "w", encoding="utf-8") as state_file:
                json.dump(merged, state_file, indent=2, sort_keys=True)
        except Exception as e:
            self.log_message(f"WARNING: Could not save processing state: {e}")

    def output_state_key(self, path):
        """Return a stable registry key for an output path."""
        if not path:
            return ""
        try:
            if self.dialog.project_folder:
                return os.path.relpath(
                    path,
                    workspace_folder(self.dialog.project_folder),
                ).replace("\\", "/")
        except ValueError:
            pass
        return os.path.abspath(path).replace("\\", "/")

    def mark_output_prepared(self, path, name=None, loaded=False, algorithm=None):
        """Record that an output file has been prepared."""
        if not path:
            return

        if not os.path.exists(path):
            return

        key = self.output_state_key(path)
        if not key:
            return

        outputs = self.processing_state.setdefault("outputs", {})
        entry = outputs.get(key, {})
        entry.update({
            "path": key,
            "absolute_path": os.path.abspath(path),
            "exists": True,
            "loaded": bool(loaded),
            "updated_at": utc_timestamp(),
        })
        if name:
            entry["name"] = name
        if algorithm:
            entry["algorithm"] = algorithm

        outputs[key] = entry
        self.save_processing_state()

    def is_output_prepared(self, path):
        """Return True when an output is recorded and present on disk."""
        if not path:
            return False

        if os.path.exists(path):
            key = self.output_state_key(path)
            if key not in self.processing_state.get("outputs", {}):
                self.mark_output_prepared(path, name=os.path.basename(path))
            return True

        key = self.output_state_key(path)
        entry = self.processing_state.get("outputs", {}).get(key)
        if entry:
            entry["exists"] = False
            self.save_processing_state()
        return False

    def record_processing_outputs(self, algorithm, params, result):
        """Record file outputs declared by a processing call."""
        if not result:
            return

        output_values = []
        for key, value in params.items():
            if key.upper().startswith("OUTPUT") and isinstance(value, str):
                output_values.append(value)

        if isinstance(result, dict):
            for key, value in result.items():
                if key.upper().startswith("OUTPUT") and isinstance(value, str):
                    output_values.append(value)

        for output_path in output_values:
            if output_path and output_path != "TEMPORARY_OUTPUT":
                self.mark_output_prepared(
                    output_path,
                    name=os.path.basename(output_path),
                    loaded=False,
                    algorithm=algorithm
                )

    def mark_workflow_status(self, workflow, status, message="", **metadata):
        """Record a project-local workflow status such as execute-all completion."""
        if not workflow:
            return

        timestamp = utc_timestamp()
        workflows = self.processing_state.setdefault("workflows", {})
        entry = workflows.get(workflow, {})
        entry.update({
            "status": status,
            "updated_at": timestamp,
        })
        if status == "running":
            entry["started_at"] = timestamp
            entry.pop("completed_at", None)
            entry.pop("failed_at", None)
        elif status == "completed":
            entry["completed_at"] = timestamp
        elif status == "failed":
            entry["failed_at"] = timestamp

        if message:
            entry["message"] = message
        elif status == "running":
            entry.pop("message", None)

        for key, value in metadata.items():
            if value is not None:
                entry[key] = value

        workflows[workflow] = entry
        self.save_processing_state()

    def save_domain_plan(self, plan):
        """Record each domain's polygon and target DEM for Morphology Setup."""
        self.processing_state["domains"] = [
            {
                "domain_id": entry["domain_id"],
                "outlet_id": entry["outlet_id"],
                "name": entry["name"],
                "polygon": entry["polygon"],
                "directory": entry["directory"],
                "dem_path": entry["dem_path"],
            }
            for entry in plan
        ]
        self.save_processing_state()
        return self.processing_state["domains"]

    def saved_domain_plan(self):
        """Return the recorded domain plan."""
        plan = self.processing_state.get("domains")
        return list(plan) if isinstance(plan, list) else []

    def workflow_status(self, workflow):
        """Return a saved workflow status entry."""
        return self.processing_state.get("workflows", {}).get(workflow, {})

    def save_grid_contract(self, l0_header, l2_header, multiplier):
        """Record the validated L0/L2 headers and multiplier for later resumes."""
        from ...grid_resolution import validate_l0_l2_alignment

        ratio = int(multiplier)
        validate_l0_l2_alignment(l0_header, l2_header, ratio)
        self.processing_state["grid"] = {
            "l0_header": dict(l0_header),
            "l2_header": dict(l2_header),
            "l2_ratio_to_l0": ratio,
            "updated_at": utc_timestamp(),
        }
        self.save_processing_state()
        return self.processing_state["grid"]

    def saved_grid_contract(self):
        """Return the saved L0/L2 grid contract, or None when it is unusable."""
        from ...grid_resolution import validate_l0_l2_alignment

        grid = self.processing_state.get("grid") or {}
        l0_header = grid.get("l0_header")
        l2_header = grid.get("l2_header")
        ratio = grid.get("l2_ratio_to_l0")
        if not isinstance(l0_header, dict) or not isinstance(l2_header, dict):
            return None
        try:
            validate_l0_l2_alignment(l0_header, l2_header, int(ratio))
        except (TypeError, ValueError) as error:
            self.log_message(
                f"WARNING: Discarding the saved L0/L2 grid contract: {error}"
            )
            return None
        return grid

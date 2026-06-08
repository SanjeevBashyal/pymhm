# -*- coding: utf-8 -*-
"""Project-local registry for prepared morphology outputs."""
from ..common import (
    os,
    json,
    processing,
)
from ...time_utils import utc_timestamp


class ProcessingStateMixin:
    """Project-local registry for prepared morphology outputs."""

    def processing_state_path(self):
        """Return the project-local processing state file path."""
        if not self.dialog.project_folder:
            return None
        return os.path.join(
            self.dialog.project_folder,
            self.processing_state_filename
        )

    def load_processing_state(self):
        """Load the processing output registry for this project."""
        state_path = self.processing_state_path()
        if not state_path or not os.path.exists(state_path):
            self.processing_state = {"version": 1, "outputs": {}}
            return

        try:
            with open(state_path, "r", encoding="utf-8") as state_file:
                state = json.load(state_file)
            if not isinstance(state, dict):
                raise ValueError("Processing state is not a JSON object.")
            state.setdefault("version", 1)
            state.setdefault("outputs", {})
            self.processing_state = state
        except Exception as e:
            self.log_message(f"WARNING: Could not read processing state: {e}")
            self.processing_state = {"version": 1, "outputs": {}}

    def save_processing_state(self):
        """Write the processing output registry to the project folder."""
        state_path = self.processing_state_path()
        if not state_path:
            return

        try:
            with open(state_path, "w", encoding="utf-8") as state_file:
                json.dump(
                    self.processing_state,
                    state_file,
                    indent=2,
                    sort_keys=True
                )
        except Exception as e:
            self.log_message(f"WARNING: Could not save processing state: {e}")

    def output_state_key(self, path):
        """Return a stable registry key for an output path."""
        if not path:
            return ""
        try:
            if self.dialog.project_folder:
                return os.path.relpath(path, self.dialog.project_folder).replace("\\", "/")
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

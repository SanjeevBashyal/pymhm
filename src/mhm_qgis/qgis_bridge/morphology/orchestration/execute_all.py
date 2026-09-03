# -*- coding: utf-8 -*-
"""Execute-all morphology workflow orchestration."""
from __future__ import annotations

from ..common import (morph_staging_folder, os,
                      project_geometry_folder)
from ....core.morphology.hydrology.outlets import configured_gauged_outlet_ids
from ..hydrology.aggregate import HydrologyMixin
from ..latlon import LatLonMixin
from ..layers import LayerProcessingMixin
from ..watershed import WatershedMixin
from ....core.executions import meteo as meteo_execution


class ExecuteAllMixin(
    WatershedMixin,
    HydrologyMixin,
    LayerProcessingMixin,
    LatLonMixin,
):
    """Execute-all morphology workflow orchestration."""

    def execute_morph_setup_processing(
            self,
            show_error_dialog=True,
            workflow_key="morph_setup") -> bool:
        """Run crop, mask, lat/lon, and ASCII export as one workflow."""
        self.log_message("\n=== Starting Morphology Setup ===")

        if not self.check_prerequisites():
            message = "Prerequisites check failed. Aborting Morphology Setup."
            self.log_message(f"ERROR: {message}")
            self.mark_workflow_status(workflow_key, "failed", message)
            return False

        self.skip_loading = True
        self.mark_workflow_status(workflow_key, "running")
        try:
            def fail(message):
                self.log_message(f"ERROR: {message}")
                self.mark_workflow_status(workflow_key, "failed", message)
                return False

            steps = meteo_execution.MORPH_SETUP_STEPS
            total = len(steps)
            for index, step in enumerate(steps, start=1):
                status_key = f"{workflow_key}_{step.key}"
                self.log_message(
                    f"\n--- Morphology Setup Step {index}/{total}: {step.label} ---")
                self.mark_workflow_status(status_key, "running")
                run = getattr(self, step.method)
                if not (run(show_error_dialog=show_error_dialog)
                        if step.takes_error_dialog else run()):
                    message = f"{step.label} failed. Aborting Morphology Setup."
                    self.mark_workflow_status(status_key, "failed", message)
                    return fail(message)
                self.mark_workflow_status(status_key, "completed")

            self.log_message("\n=== Morphology Setup Completed Successfully ===")
            self.mark_workflow_status(
                workflow_key,
                "completed",
                "Morphology Setup completed successfully.",
            )
            return True

        except Exception as e:
            message = f"Morphology Setup failed with exception: {str(e)}"
            self.log_message(f"\nERROR: {message}")
            import traceback

            self.log_message(f"Traceback: {traceback.format_exc()}")
            self.mark_workflow_status(workflow_key, "failed", message)
            if show_error_dialog:
                self.error("Error", f"Morphology Setup failed:\n{str(e)}"
                )
            return False
        finally:
            self.skip_loading = False
            self.log_message(
                "Morphology Setup finished. Prepared outputs were recorded in the project processing state."
            )

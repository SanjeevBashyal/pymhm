# -*- coding: utf-8 -*-
"""Project-local registry for prepared morphology outputs."""
from ..common import (
    os,
    json,
)
from ....core.handlers.store.paths import workspace_folder
from ....core.handlers.store.registry import available, key_for, register
from ....core.utils.time import utc_timestamp
from ....core.handlers.state import morphology as morphology_state
from ....core.session import MorphologySession


class ProcessingStateMixin:
    """Project-local registry for prepared morphology outputs."""

    def _state(self) -> MorphologySession:
        """Return the run state, creating it for objects built without __init__.

        The tests exercise a single mixin by subclassing it directly, so the
        processor's __init__ has not run and there is no session yet.
        """
        session = self.__dict__.get("_session")
        if session is None:
            # Wire the callbacks too: a session created here must be as usable
            # as one built in __init__, or its messages go nowhere.
            session = MorphologySession(
                project_folder=getattr(
                    getattr(self, "dialog", None), "project_folder", "") or "",
                log=getattr(self, "log_message", None),
                load=getattr(self, "load_layer", None),
                run=getattr(self, "run_processing_algorithm", None))
            self.__dict__["_session"] = session
        return session

    @property
    def session(self) -> MorphologySession:
        """The run state, with the project folder refreshed from the dialog."""
        session = self._state()
        # Refresh the callbacks every time rather than capturing them once: an
        # object may acquire log_message or load_layer after the session was
        # first created, and a stale None sends its messages nowhere.
        for field, attr in (("log", "log_message"), ("load", "load_layer"),
                            ("run", "run_processing_algorithm")):
            found = getattr(self, attr, None)
            if found is not None:
                setattr(session, field, found)
        dialog = getattr(self, "dialog", None)
        if dialog is None:
            return session
        session.project_folder = getattr(dialog, "project_folder", "") or ""
        # Resolve the widget-backed values once, here, so no step below has to
        # touch a widget: the CRS is a value and the selections are sources.
        get_crs = getattr(dialog, "get_crs", None)
        if get_crs is not None:
            try:
                session.crs = get_crs()
            except Exception:
                pass
        # Ask the dialog which input kinds exist rather than naming them here:
        # the spec kinds are not the spec keys ("lc" selects "land_cover"), and
        # a hard-coded list silently resolves nothing for the ones it misses.
        kinds = set(getattr(dialog, "_input_adapters", {}) or ())
        kinds.update(("dem", "pour_points"))
        for kind in kinds:
            found = self._selected_source(kind)
            if found:
                session.input_sources[kind] = found
        session.dem_source = session.source_for("dem") or session.dem_source
        session.pour_points_source = (
            session.source_for("pour_points") or session.pour_points_source)
        for name, attr in (("grid_level_headers", "grid_headers"),
                           ("current_l0_resolution", "l0_resolution"),
                           ("filled_dem_resolution_info", "dem_resolution_info")):
            getter = getattr(dialog, name, None)
            if getter is not None:
                try:
                    setattr(session, attr, getter())
                except Exception:
                    pass
        return session

    def _selected_source(self, kind: str) -> str:
        """Return the source of a selected input layer, or ''."""
        combo = getattr(self.dialog, "input_combo", None)
        if combo is None:
            return ""
        try:
            box = combo(kind)
        except Exception:
            return ""
        if box is None:
            return ""
        # A folder-backed adapter knows its own path; a QGIS layer does not.
        if hasattr(box, "source_path"):
            found = box.source_path() or ""
            if found:
                return str(found)
        layer = box.currentLayer()
        return str(layer.source() or "") if layer is not None else ""

    def processing_state_path(self):
        return morphology_state.state_path(self.session)

    def load_processing_state(self):
        return morphology_state.load(self.session)

    def save_processing_state(self):
        return morphology_state.save(self.session)

    def output_state_key(self, path):
        return morphology_state.output_key(self.session, path)

    def mark_output_prepared(self, path, name=None, loaded=False, algorithm=None):
        return morphology_state.mark_prepared(self.session, path, name=name, loaded=loaded, algorithm=algorithm)

    def is_output_prepared(self, path):
        return morphology_state.is_prepared(self.session, path)

    def record_processing_outputs(self, algorithm, params, result):
        return morphology_state.record_outputs(self.session, algorithm, params, result)

    def mark_workflow_status(self, workflow, status, message="", **metadata):
        return morphology_state.mark_workflow(self.session, workflow, status, message, **metadata)

    def workflow_status(self, workflow):
        return morphology_state.workflow_status(self.session, workflow)

    def save_domain_plan(self, plan):
        return morphology_state.save_domain_plan(self.session, plan)

    def saved_domain_plan(self):
        return morphology_state.saved_domain_plan(self.session)

    def save_grid_contract(self, l0_header, l2_header, multiplier):
        return morphology_state.save_grid_contract(self.session, l0_header, l2_header, multiplier)

    def saved_grid_contract(self):
        return morphology_state.saved_grid_contract(self.session)

    @property
    def processing_state(self):
        return self._state().processing_state

    @processing_state.setter
    def processing_state(self, value):
        self._state().processing_state = value

    @property
    def skip_loading(self):
        return self._state().skip_loading

    @skip_loading.setter
    def skip_loading(self, value):
        self._state().skip_loading = value

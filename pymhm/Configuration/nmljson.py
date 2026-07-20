# -*- coding: utf-8 -*-
"""Project-local namelist JSON file profiles."""
from __future__ import annotations

import json
import os
import tempfile
from typing import Any

from ..time_utils import utc_timestamp
from .paths import output_filename


NMLJSON_FILENAME = "nmljson.json"
PROFILE_KINDS = ("mhm", "parameters", "outputs")


def atomic_write_text(path: str, content: str) -> None:
    """Replace a UTF-8 text file only after its complete content is written."""
    directory = os.path.dirname(path) or "."
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".nmljson-", dir=directory)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


class NmlJsonStore:
    """Read and write the three converter-ready namelist file profiles."""

    def __init__(self, dialog: Any) -> None:
        self.dialog = dialog

    def path(self) -> str | None:
        if not self.dialog.project_folder:
            return None
        return os.path.join(self.dialog.project_folder, NMLJSON_FILENAME)

    def empty(self, mhm_version: str = "") -> dict[str, Any]:
        return {
            "format_version": 1,
            "mhm_version": mhm_version,
            "file_profiles": {
                kind: self._empty_profile(kind, mhm_version)
                for kind in PROFILE_KINDS
            },
        }

    def _empty_profile(self, kind: str, mhm_version: str) -> dict[str, Any]:
        return {
            "format_version": 1,
            "profile": kind,
            "output_file": output_filename(kind, mhm_version),
            "namelists": [],
            "values": {},
        }

    def load(self, mhm_version: str = "") -> dict[str, Any]:
        path = self.path()
        if not path or not os.path.exists(path):
            return self.empty(mhm_version)
        try:
            with open(path, "r", encoding="utf-8") as handle:
                document = json.load(handle)
            if not isinstance(document, dict):
                raise ValueError("nmljson root is not an object")
            version = str(document.get("mhm_version") or mhm_version)
            document.setdefault("format_version", 1)
            document["mhm_version"] = version
            profiles = document.setdefault("file_profiles", {})
            if not isinstance(profiles, dict):
                raise ValueError("nmljson file_profiles is not an object")
            for kind in PROFILE_KINDS:
                profile = profiles.setdefault(
                    kind, self._empty_profile(kind, version))
                profile.setdefault("profile", kind)
                profile.setdefault("format_version", 1)
                profile.setdefault("output_file", output_filename(kind, version))
                profile.setdefault("namelists", [])
                profile.setdefault("values", {})
            return document
        except Exception as exc:
            self.dialog.log_message(f"WARNING: Could not read nmljson: {exc}")
            return self.empty(mhm_version)

    def save(self, document: dict[str, Any]) -> None:
        path = self.path()
        if not path:
            return
        document["format_version"] = 1
        document["updated_at"] = utc_timestamp()
        atomic_write_text(
            path,
            json.dumps(document, indent=2, sort_keys=False) + "\n",
        )

    def profile(self, document: dict[str, Any], kind: str) -> dict[str, Any]:
        return document.setdefault("file_profiles", {}).setdefault(
            kind,
            self._empty_profile(kind, str(document.get("mhm_version") or "")),
        )

    def set_profile(
            self,
            document: dict[str, Any],
            kind: str,
            values: dict[str, Any],
            mhm_version: str,
            editor_values: dict[str, Any] | None = None,
            dimensions: dict[str, int] | None = None) -> dict[str, Any]:
        document["mhm_version"] = mhm_version
        profile = self.profile(document, kind)
        profile.update({
            "format_version": 1,
            "profile": kind,
            "output_file": output_filename(kind, mhm_version),
            "namelists": list(values),
            "values": values,
        })
        if dimensions:
            profile["dimensions"] = dimensions
        else:
            profile.pop("dimensions", None)
        if editor_values is not None:
            profile["editor_values"] = editor_values
        return document

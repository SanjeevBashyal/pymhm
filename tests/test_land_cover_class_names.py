"""Tests for land-cover labels used by elevation-band reports."""

import pandas as pd

from mhm_qgis import standalone

standalone.install(force=True)

from mhm_qgis.Morphology.layers.land_cover_class_names import (  # noqa: E402
    LandCoverClassNameMixin,
)


class _Dialog:
    def categorical_lookup_config(self, _kind):
        return {
            "lookup_table": "lookup.csv",
            "mapping_field": "source_value",
            "class_field": "class_id",
        }


class _Reader(LandCoverClassNameMixin):
    dialog = _Dialog()

    def log_message(self, _message):
        pass


def test_names_are_keyed_by_formatted_class_field(monkeypatch):
    # Patch the lookup API, which is the boundary the code under test crosses.
    # Reaching into mhm_tools.common directly coupled this to where mhm-tools
    # happened to keep read_lookup_table, and it moved.
    from mhm_qgis.core.handlers import lookup

    table = pd.DataFrame(
        {
            "source_value": [10, 20],
            "class_id": [1, 2],
            "name": ["Forest", "Urban"],
        }
    )
    monkeypatch.setattr(
        lookup, "read", lambda _path: table
    )

    assert _Reader()._read_land_cover_class_names() == {
        1: "Forest",
        2: "Urban",
    }

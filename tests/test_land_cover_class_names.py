"""Tests for land-cover labels used by elevation-band reports."""

import pandas as pd

from mhm_qgis.Morphology.layers.land_cover_class_names import \
    LandCoverClassNameMixin


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
    from mhm_tools.common import format_data

    table = pd.DataFrame(
        {
            "source_value": [10, 20],
            "class_id": [1, 2],
            "name": ["Forest", "Urban"],
        }
    )
    monkeypatch.setattr(format_data, "read_lookup_table", lambda _path: table)

    assert _Reader()._read_land_cover_class_names() == {
        1: "Forest",
        2: "Urban",
    }

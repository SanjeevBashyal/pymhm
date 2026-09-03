"""Optional land-cover class names used by elevation-band reports."""

from __future__ import annotations

from ....core.handlers import lookup


class LandCoverClassNameMixin:
    """Read reporting labels from the active land-cover lookup table."""

    def _read_land_cover_class_names(self):
        config = self.dialog.categorical_lookup_config("lc")
        if not config:
            return {}

        try:
            table = lookup.read(config["lookup_table"])
            fields = {
                lookup.normalize_key(column): column for column in table.columns
            }
            class_field = fields[lookup.normalize_key(config["class_field"])]
            name = next(
                (
                    fields[field]
                    for field in (
                        "classname",
                        "name",
                        "label",
                        "description",
                        "landcovername",
                        "type",
                    )
                    if field in fields
                ),
                None,
            )
            if name is None:
                return {}
            return {
                int(float(row[class_field])): str(row[name]).strip()
                for _, row in table.iterrows()
                if str(row[name]).strip()
            }
        except Exception as error:
            self.log_message(
                f"WARNING: Could not read land-cover class names: {error}"
            )
            return {}


__all__ = ["LandCoverClassNameMixin"]

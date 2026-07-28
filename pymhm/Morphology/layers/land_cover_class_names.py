"""Optional land-cover class names used by elevation-band reports."""

from __future__ import annotations


def _normalized(value):
    text = str(value).strip().lstrip("*").split("[", 1)[0]
    return "".join(character.lower() for character in text if character.isalnum())


class LandCoverClassNameMixin:
    """Read reporting labels from the active land-cover lookup table."""

    def _read_land_cover_class_names(self):
        config = self.dialog.categorical_lookup_config("lc")
        if not config:
            return {}

        from ...mhm_tools_adapter import read_categorical_lookup_table

        try:
            table = read_categorical_lookup_table(config["lookup_table"])
            fields = {_normalized(column): column for column in table.columns}
            class_field = fields[_normalized(config["class_field"])]
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

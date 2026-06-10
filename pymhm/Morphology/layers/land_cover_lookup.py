# -*- coding: utf-8 -*-
"""Compatibility aggregate for land-cover lookup helpers."""
from .lookup_fields import LookupFieldMixin
from .land_cover_lookup_table import LandCoverLookupTableMixin
from .land_cover_class_names import LandCoverClassNameMixin


class LandCoverLookupMixin(
    LandCoverLookupTableMixin,
    LandCoverClassNameMixin,
    LookupFieldMixin,
):
    """Compatibility aggregate for land-cover lookup helpers."""

    pass

"""Tools to prepare, process, evaluate, and visualize mHM data.

.. toctree::
   :hidden:

   self

Subpackages
===========

Built-in processing and tool functions.

.. autosummary::
   :toctree: api
   :caption: Subpackages

   common
   data_processing
   evaluation
   legacy_tools
   mhm_v5_v6_converter
   setup_creation
   utility_tools
   visualisation
"""

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "not_available"

from . import (
    common,
    data_processing,
    evaluation,
    legacy_tools,
    mhm_v5_v6_converter,
    setup_creation,
    utility_tools,
    visualisation,
)

__all__ = ["__version__"]
__all__ += [
    "common",
    "data_processing",
    "evaluation",
    "legacy_tools",
    "mhm_v5_v6_converter",
    "setup_creation",
    "utility_tools",
    "visualisation",
]

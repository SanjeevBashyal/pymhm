"""General utility tools for mHM workflows.

Submodules and public entry points are imported lazily following the package
pattern used by the legacy ``pre`` and ``post`` packages.
"""

import importlib

_ATTR_EXPORTS = {
    "link_folder_tree": ("link_folder_tree", "link_folder_tree"),
}

__all__ = [
    "link_folder_tree",
]


def __getattr__(name):
    if name in _ATTR_EXPORTS:
        module_name, attr_name = _ATTR_EXPORTS[name]
        module = importlib.import_module(f".{module_name}", __name__)
        attr = getattr(module, attr_name)
        globals()[name] = attr
        return attr
    error_msg = f"module '{__name__}' has no attribute '{name}'"
    raise AttributeError(error_msg)


def __dir__():
    return sorted(set(globals().keys()) | set(__all__))

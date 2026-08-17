from importlib import import_module
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("senid")
except PackageNotFoundError:
    __version__ = "0+unknown"


_MODULES = {
    # Short public aliases
    "pp": "preprocessing",
    "it": "intrinsic",
    "et": "extrinsic",
    "tl": "tools",
    "pl": "plotting",

    # Full module names
    "preprocessing": "preprocessing",
    "intrinsic": "intrinsic",
    "extrinsic": "extrinsic",
    "tools": "tools",
    "plotting": "plotting",
    "spatial": "spatial",
}


__all__ = [
    "__version__",
    *_MODULES.keys(),
]


def __getattr__(name: str):
    module_name = _MODULES.get(name)

    if module_name is not None:
        module = import_module(f"senid.{module_name}")

        # Cache it so __getattr__ is only invoked once.
        globals()[name] = module

        return module

    raise AttributeError(
        f"module 'senid' has no attribute {name!r}"
    )
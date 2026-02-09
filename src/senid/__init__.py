"""
SenID: identify cell-intrinsic and cell-extrinsic features of senescence and communication.
"""

from __future__ import annotations
from importlib import import_module

__all__ = [
    "__version__",
]

# Keep this manual or wire it later via importlib.metadata
__version__ = "0.1.1"

def __getattr__(name: str):
    if name in {"intrinsic", "extrinsic", "preprocessing", "plotting", "tools", "spatial"}:
        return import_module(f"senid.{name}")
    raise AttributeError(f"module 'senid' has no attribute {name!r}")
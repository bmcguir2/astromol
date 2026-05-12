"""Public package metadata for astromol."""

from __future__ import annotations

from importlib import metadata

try:
    __version__ = metadata.version("astromol")
except metadata.PackageNotFoundError:  # pragma: no cover - source checkout fallback
    __version__ = "0.1.0.dev0"

__all__ = ["__version__"]

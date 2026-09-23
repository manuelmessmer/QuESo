"""Python interface for QuESo embedded finite-element model generation."""

import sys

MIN_PYTHON_VERSION = (3, 9)
if sys.version_info < MIN_PYTHON_VERSION:
    version = ".".join(map(str, MIN_PYTHON_VERSION))
    raise RuntimeError(f"pyqueso requires Python {version} or higher.")

try:
    from . import _core
except ImportError as exc:
    raise ImportError(f"Failed to import pyqueso._core: {exc}") from exc

from .model import Model

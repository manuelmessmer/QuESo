import sys
from pathlib import Path

# Enforce minimum Python version
MIN_PYTHON_VERSION = (3, 9)
if sys.version_info < MIN_PYTHON_VERSION:
    raise RuntimeError(f"QuESo requires Python {'.'.join(map(str, MIN_PYTHON_VERSION))} or higher.")

# The module is located in the same directory as the present file.
QUESO_PYTHON_LIBS_PATH = Path(__file__).resolve().parent
if QUESO_PYTHON_LIBS_PATH not in sys.path:
    sys.path.append(QUESO_PYTHON_LIBS_PATH)

# Import _QuESoPythonModule modules
try:
    from ._QuESoPythonModule import PrintLogo
    PrintLogo()
    from ._QuESoPythonModule import *
except ImportError as e:
    raise ImportError(f"Failed to import _QuESoPythonModule {e}")



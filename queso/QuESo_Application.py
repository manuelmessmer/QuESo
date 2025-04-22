import os
import sys

# Enforce minimum Python version
MIN_PYTHON = (3, 5)
if sys.version_info < MIN_PYTHON:
    raise RuntimeError(f"QuESo requires Python {'.'.join(map(str, MIN_PYTHON))} or higher.")

class Paths:
    # Base install path is the parent directory of this script's location
    queso_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    queso_libs = os.path.join(queso_install_path, "libs")

# Append QuESo library path to sys.path if not already included
if Paths.queso_libs not in sys.path:
    sys.path.append(Paths.queso_libs)

# Import QuESo application modules
try:
    from QuESo_Application import PrintLogo
    PrintLogo()
    from QuESo_Application import *
except ImportError as e:
    raise ImportError(f"Failed to import QuESo_Application modules: {e}")



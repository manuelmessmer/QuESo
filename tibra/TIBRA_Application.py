# Application dependent names and paths
import os
import sys

if sys.version_info < (3, 5):
    raise Exception("Kratos only supports Python version 3.5 and above")

class Paths(object):
    tibra_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    tribera_libs = os.path.join(tibra_install_path,"libs")

# import core library (Kratos.so)
sys.path.append(Paths.tribera_libs)
from TIBRA_Application import *



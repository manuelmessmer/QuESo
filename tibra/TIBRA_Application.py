# Application dependent names and paths
import os
import sys

if sys.version_info < (3, 5):
    raise Exception("TIBRA only supports Python version 3.5 and above")

class Paths(object):
    tibra_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    tibra_libs = os.path.join(tibra_install_path,"libs")

# import core library (Kratos.so)
sys.path.append(Paths.tibra_libs)
from TIBRA_Application import *



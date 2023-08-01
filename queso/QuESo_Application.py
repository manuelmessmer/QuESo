# Application dependent names and paths
import os
import sys

if sys.version_info < (3, 5):
    raise Exception("QuESo only supports Python version 3.5 and above")

class Paths(object):
    queso_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    queso_libs = os.path.join(queso_install_path,"libs")

# import core library (Kratos.so)
sys.path.append(Paths.queso_libs)
from QuESo_Application import *



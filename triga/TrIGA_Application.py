# Application dependent names and paths
import os
import sys

if sys.version_info < (3, 5):
    raise Exception("Kratos only supports Python version 3.5 and above")

class Paths(object):
    cgal_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    print(cgal_install_path)
    cgal_libs = os.path.join(cgal_install_path,"libs")
    print(cgal_libs)

# import core library (Kratos.so)
sys.path.append(Paths.cgal_libs)
from TrIGA_Application import *



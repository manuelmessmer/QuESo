# Application dependent names and paths
import os
import sys

if sys.version_info < (3, 5):
    raise Exception("TrIGA only supports Python version 3.5 and above")

class Paths(object):
    TrIGA_install_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    print(TrIGA_install_path)
    TrIGA_libs = os.path.join(TrIGA_install_path,"libs")
    print(TrIGA_libs)

# import core library (TrIGA.so)
sys.path.append(Paths.TrIGA_libs)
from TrIGA_Application import *



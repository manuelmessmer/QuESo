# Trivariate Isogeometric B-Rep Analysis - TIBRA 

[![Build Status](https://github.com/manuelmessmer/TIBRA/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/manuelmessmer/TIBRA/actions)

TIBRA is an IGA preprocessor to generate analysis-ready models from volumetric B-Reps. The B-Rep is embedded into a trivariate B-Spline cuboid defined by a bounding box (`lower_point`, `upper_point`). The user can specify the polynomial degree (`p=1..4`) and the number of knot spans in each spatial direction. TIBRA reads an STL file and computes the integration points required for an FE analysis. An interface to the FE framework [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) is provided. More information about the theoretical background of TIBRA can be found in [1]. Please do not hesitate to contact me if you have any questions about TIBRA.

This project is still under development. TIBRA is written in C++ and has a user-friendly Python interface.

Before using TIBRA, you might want to check if all unit tests pass: `./run_tests all`

Different examples how to use TIBRA can be found in [Examples](https://github.com/manuelmessmer/TIBRA/tree/main/examples).

Input: Solid B-Rep Model             |  Output: IGA Model ready for FE Analysis
:-------------------------:|:-------------------------:
![](https://github.com/manuelmessmer/TIBRA/blob/main/docs/brep.png)  |  ![](https://github.com/manuelmessmer/TIBRA/blob/main/docs/iga_model.png)

## Installation (Linux)
Required packages:
```console
sudo apt-get install python3-dev gcc g++ cmake libboost-test-dev
```
Required Python modules:
```console
sudo pip3 install numpy scipy
```
Installation ([configure.sh](https://github.com/manuelmessmer/TIBRA/tree/main/configure.sh)):
```console
sh configure.sh
```
To run TIBRA make sure LD_LIBRARY_PATH and PYTHONPATH are pointing to your TIBRA folder:
```console
export PYTHONPATH=/path/to/TIBRA/
```
```console
export LD_LIBRARY_PATH=/path/to/TIBRA/libs/
```

## How to cite TIBRA?
Please use the following references when citing TIBRA in your work.
- [1] Manuel Meßmer, Tobias Teschemacher, Lukas F. Leidinger, Roland Wüchner, Kai-Uwe Bletzinger, Efficient CAD-integrated isogeometric analysis of trimmed solids, Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
- [2] Manuel Meßmer, Lukas F. Leidiner, Stefan Hartmann, ..., Kai-Uwe Bletzinger, Isogeometric Analysis on Trimmed Solids: A B-Spline-Based Approach Focusing on Explicit Dynamics, 13th European LS-DYNA Conference, Ulm, Germany, 2021. [Meßmer et al. 2022](https://www.researchgate.net/publication/357053531_Isogeometric_Analysis_on_Trimmed_Solids_A_B-Spline-Based_Approach_Focusing_on_Explicit_Dynamics).
- [3] Manuel Meßmer, TIBRA, https://github.com/manuelmessmer/TIBRA.


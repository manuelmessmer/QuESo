# Trivariate Isogeometric B-Rep Analysis - TIBRA 

TIBRA is an IGA preprocessor to generate analysis suitable models from volumetric B-Reps. The B-Rep is embedded into a trivariate B-Spline cuboid defined by a bounding box (`lower_point`, `upper_point`). The user can specify the polynomial degree (`p=1..4`) and the number of knot spans in each spatial direction. TIBRA reads an STL file and computes the integration points required for a FE analysis. An interface to the FE framework Kratos Multiphysics is provided. More information about the theoretical background of TIBRA can be found in [1]. Please do not hesitate to contact me if you have any questions about TIBRA.

This project is still under development.

Before using TIBRA, please check if all unit tests pass: `./run_all_test cpp`

The python tests require a Kratos Multiphysics installation: `./run_all_test py`

A simple example how to use TIBRA can be found here: [Example](https://github.com/manuelmessmer/TIBRA/tree/main/examples/cantilever).

## Installation
Required dependencies (C++):
- CMAKE (minimum 3.15)
- BOOST Unit Test Framework
- OpenMP
- CGAL - https://www.cgal.org/ (minimum version 5.3.2) 

Required Python modules:
- json

Execute `sh configure.sh` to install TIBRA.

## How to cite TIBRA?
Please use the following references when citing Kratos in your work.
- [1] Manuel Meßmer, Tobias Teschemacher, Lukas F. Leidinger, Roland Wüchner, Kai-Uwe Bletzinger, Efficient CAD-integrated isogeometric analysis of trimmed solids, Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
- [2] Manuel Meßmer, Lukas F. Leidiner, Stefan Hartmann, ..., Kai-Uwe Bletzinger, Isogeometric Analysis on Trimmed Solids: A B-Spline-Based Approach Focusing on Explicit Dynamics, 13th European LS-DYNA Conference, Ulm, Germany, 2021. [Meßmer et al. 2022](https://www.researchgate.net/publication/357053531_Isogeometric_Analysis_on_Trimmed_Solids_A_B-Spline-Based_Approach_Focusing_on_Explicit_Dynamics).
- [3] Manuel Meßmer, TIBRA, https://github.com/manuelmessmer/TIBRA.


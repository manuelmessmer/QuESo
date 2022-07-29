# Trivariate/Trimmed Isogeometric Analysis - TrIGA 

TrIGA is an IGA preprocessor to generate analysis suitable models from volumetric B-Reps. The B-Rep is embedded into a trivariate B-Spline cuboid defined by a bounding box (`lower_point`, `upper_point`). The user can specify the polynomial degree (`p=1..4`) and the number of knot spans in each spatial direction. TrIGA reads an STL file and computes the integration points required for a FE analysis. An interface to the FE framework Kratos Multiphysics is provided. 

This project is still under development.

Before using TrIGA, please check if all unit tests pass: `./run_all_test cpp`

The python tests require a Kratos Multiphysics installation: `./run_all_test py`

## Installation
Required dependencies (C++):
- CMAKE (minimum 3.15)
- BOOST Unit Test Framework
- OpenMP
- CGAL - https://www.cgal.org/ (minimum version 5.3.2) 

Required Python modules:
- json

Execute `sh configure.sh` to install TrIGA.

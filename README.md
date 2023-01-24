# Trivariate Isogeometric B-Rep Analysis - TIBRA 

[![License: MIT](https://img.shields.io/badge/License-BSD4-green.svg)](https://github.com/manuelmessmer/TIBRA/blob/main/LICENSE) [![C++][c++-image]][c++standard] [![Build Status](https://github.com/manuelmessmer/TIBRA/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/manuelmessmer/TIBRA/actions)

[c++-image]: https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B
[c++standard]: https://isocpp.org/std/the-standard

TIBRA is an IGA preprocessor to generate analysis-ready models from volumetric B-Reps. The B-Rep is embedded into a trivariate B-Spline cuboid defined by a bounding box. The user can specify the polynomial degree and the number of knot spans in each spatial direction. TIBRA reads an STL file and computes the integration points, control points and knot vectors required for an FE analysis. An interface to the FE framework [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) is provided. TIBRA is written in C++ and has a user-friendly Python interface.
More information about the theoretical background of TIBRA can be found in [1]. Please do not hesitate to contact me if you have any questions about TIBRA.

If you want to learn how to use TIBRA, check out the [Wiki](https://github.com/manuelmessmer/TIBRA/wiki/Getting-Started). Additionally, there are several examples in [Examples](https://github.com/manuelmessmer/TIBRA/tree/main/examples).

Input: Solid B-Rep Model             |  Output: IGA Model ready for FE Analysis
:-------------------------:|:-------------------------:
![](https://github.com/manuelmessmer/TIBRA/blob/main/docs/brep.png)  |  ![](https://github.com/manuelmessmer/TIBRA/blob/main/docs/iga_model.png)

## How to cite TIBRA?
Please use the following references when citing TIBRA in your work.
- [1] Manuel Meßmer, Tobias Teschemacher, Lukas F. Leidinger, Roland Wüchner, Kai-Uwe Bletzinger, Efficient CAD-integrated isogeometric analysis of trimmed solids, Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
- [2] Manuel Meßmer, Lukas F. Leidiner, Stefan Hartmann, ..., Kai-Uwe Bletzinger, Isogeometric Analysis on Trimmed Solids: A B-Spline-Based Approach Focusing on Explicit Dynamics, 13th European LS-DYNA Conference, Ulm, Germany, 2021. [Meßmer et al. 2022](https://www.researchgate.net/publication/357053531_Isogeometric_Analysis_on_Trimmed_Solids_A_B-Spline-Based_Approach_Focusing_on_Explicit_Dynamics).
- [3] Manuel Meßmer, TIBRA, https://github.com/manuelmessmer/TIBRA.


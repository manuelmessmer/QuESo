# QuESo - Quadrature for Embedded Solids

[![License: MIT](https://img.shields.io/badge/License-BSD4-green.svg)](https://github.com/manuelmessmer/QuESo/blob/main/LICENSE) [![C++][c++-image]][c++standard] [![Build Status](https://github.com/manuelmessmer/QuESo/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/manuelmessmer/QuESo/actions)

[c++-image]: https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B
[c++standard]: https://isocpp.org/std/the-standard

:exclamation::exclamation::exclamation: **I am currently in the process of renaming this repository. You might find some info that needs to be updated.** :exclamation::exclamation::exclamation:

QuESo is an embedding preprocessor to generate analysis-ready FE models from arbitrarily complex B-Reps, given as an STL file. The STL is embedded into the background mesh defined by a simple bounding box. The user can specify the polynomial degree and the number of elements in each spatial direction. QuESo reads the STL file and computes the integration points required for an FE analysis. Both classical $C^0$ finite elements and isogeometric elements defined on trivariate B-Spline domains are supported. An interface to the FE framework [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) is provided. QuESo is written in C++ and has a user-friendly Python interface. More information about the theoretical background of QuESo can be found in [1]. Please do not hesitate to contact me with questions about QuESo.

If you want to learn how to use QuESo, check out the [Wiki](https://github.com/manuelmessmer/QuESo/wiki/Getting-Started). Additionally, there are several examples in [Examples](https://github.com/manuelmessmer/QuESo/tree/main/examples).

Input: Solid B-Rep Model (STL)             |  Output: Quadrature Rules
:-------------------------:|:-------------------------:
![](https://github.com/manuelmessmer/QuESo/blob/main/docs/brep.png)  |  ![](https://github.com/manuelmessmer/QuESo/blob/main/docs/iga_model.png)

## How to cite QuESo?
Please use the following references when citing QuESo in your work.
- [1] Manuel Meßmer, Tobias Teschemacher, Lukas F. Leidinger, Roland Wüchner, Kai-Uwe Bletzinger, Efficient CAD-integrated isogeometric analysis of trimmed solids, Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
- [2] Manuel Meßmer, Lukas F. Leidiner, Stefan Hartmann, ..., Kai-Uwe Bletzinger, Isogeometric Analysis on Trimmed Solids: A B-Spline-Based Approach Focusing on Explicit Dynamics, 13th European LS-DYNA Conference, Ulm, Germany, 2021. [Meßmer et al. 2022](https://www.researchgate.net/publication/357053531_Isogeometric_Analysis_on_Trimmed_Solids_A_B-Spline-Based_Approach_Focusing_on_Explicit_Dynamics).
- [3] Manuel Meßmer, QuESo, https://github.com/manuelmessmer/QuESo.


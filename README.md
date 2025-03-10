<p align=center><img height="72.125%" width="72.125%" src="https://github.com/manuelmessmer/QuESo/blob/main/docs/logo.png"></p>

[![License: BSD4](https://img.shields.io/badge/License-BSD4-yellow.svg)](https://github.com/manuelmessmer/QuESo/blob/main/LICENSE) [![C++][c++-image]][c++standard] 
[![CI](https://github.com/manuelmessmer/QuESo/actions/workflows/ci.yml/badge.svg)](https://github.com/manuelmessmer/QuESo/actions)
[![Docker](https://img.shields.io/badge/Docker-Repository-blue?style=flat-square)](https://github.com/manuelmessmer/docker-queso)
![GitHub stars](https://img.shields.io/github/stars/manuelmessmer/QuESo?style=social)

[c++-image]: https://img.shields.io/badge/C++-17-blue.svg?style=flat&logo=c%2B%2B
[c++standard]: https://isocpp.org/std/the-standard

QuESo (formerly TIBRA) is a preprocessor to generate analysis-ready embedded finite element models from arbitrarily complex solids described in boundary representation (B-Rep). The geometry is processed as a standard STL file. QuESo is designed to provide highly efficient quadrature rules that can be used in any embedded or immersed boundary method. To this end, the STL model is embedded into a background mesh defined by a regular grid, i.e., with hexahedral integration domains (finite elements). The user can specify the polynomial degree and the number of elements in each spatial direction. QuESo reads the STL file and computes the integration points required for a subsequent FE analysis. Both classical $C^0$ finite elements and isogeometric elements defined on trivariate B-Spline domains are supported. An interface to the FE framework [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) is provided. The integration points have the following characteristics:
* All integration weights are positive.
* Point locations are restricted to the material domain.
* The number of points per cut element is always $n \leq (p+1)^3$.
  
QuESo is written in C++ and has a user-friendly Python interface. If you want to learn how to use QuESo, check out the [**Wiki**](https://github.com/manuelmessmer/QuESo/wiki). Additionally, there are several examples in [Examples](https://github.com/manuelmessmer/QuESo/tree/main/examples). Please do not hesitate to contact me with questions about QuESo.

![](https://github.com/manuelmessmer/QuESo/blob/main/docs/input_output.png) 

## Special Thanks To
* Lester Hedges for the [AABB tree](https://github.com/lohedges/aabbcc)
* Mike Lapshin for the [NNLS solver](https://github.com/mlapshin/nnls)
* [pybind11](https://github.com/pybind/pybind11) for exposing C++ to python
* [Boost](https://www.boost.org/users/download/) for the C++ unit test framework
  
## How to cite QuESo?
Please use the following references when citing QuESo in your work.

* Manuel Meßmer, Tobias Teschemacher, Lukas F. Leidinger, Roland Wüchner, Kai-Uwe Bletzinger, Efficient CAD-integrated isogeometric analysis of trimmed solids, Comput. Methods Appl. Mech. Engrg.  400 (2022) 115584, [https://doi.org/10.1016/j.cma.2022.115584](https://doi.org/10.1016/j.cma.2022.115584).

* Manuel Meßmer, Stefan Kollmannsberger, Roland Wüchner, Kai-Uwe Bletzinger, Robust numerical integration of embedded solids described in boundary representation, Comput. Methods Appl. Mech. Engrg.  419 (2024) 116670, [https://doi.org/10.1016/j.cma.2023.116670](https://doi.org/10.1016/j.cma.2023.116670).
  
* Manuel Meßmer, Lukas F. Leidinger, Stefan Hartmann, ..., Kai-Uwe Bletzinger, Isogeometric Analysis on Trimmed Solids: A B-Spline-Based Approach Focusing on Explicit Dynamics, 13th European LS-DYNA Conference, Ulm, Germany, 2021. [Meßmer et al. 2022](https://www.researchgate.net/publication/357053531_Isogeometric_Analysis_on_Trimmed_Solids_A_B-Spline-Based_Approach_Focusing_on_Explicit_Dynamics).

* Manuel Meßmer, QuESo, https://github.com/manuelmessmer/QuESo.

QuESo has also been used in the scope of shape optimization of complex embedded solid structures.


* Manuel Meßmer, Reza Najian Asl, Stefan Kollmannsberger, Roland Wüchner, Kai-Uwe Bletzinger, Shape optimization of embedded solids using implicit Vertex-Morphing, Comput. Methods Appl. Mech. Engrg.  426 (2024) 116999, [https://doi.org/10.1016/j.cma.2024.116999](https://doi.org/10.1016/j.cma.2024.116999)
  
It would be great to hear from you if QuESo was of use in your project.




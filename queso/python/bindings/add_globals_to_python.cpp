//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

/// Project inlcudes
#include "queso/python/bindings/add_globals_to_python.h"
// To export
#include "queso/includes/define.hpp"

namespace queso::python {

namespace py = pybind11;

void AddGlobalsToPython(pybind11::module& m)
{

    /// Export enum IntegrationMethod
    py::enum_<IntegrationMethod>(m, "IntegrationMethod")
        .value("GAUSS", IntegrationMethod::gauss, "Full Gauss quadrature.")
        .value("GAUSS_REDUCED_1", IntegrationMethod::gauss_reduced_1, "First reduced Gauss rule.")
        .value("GAUSS_REDUCED_2", IntegrationMethod::gauss_reduced_2, "Second reduced Gauss rule.")
        .value("GGQ_OPTIMAL", IntegrationMethod::ggq_optimal, "Optimal generalized Gaussian quadrature.")
        .value("GGQ_REDUCED_1", IntegrationMethod::ggq_reduced_1, "First reduced GGQ rule.")
        .value("GGQ_REDUCED_2", IntegrationMethod::ggq_reduced_2, "Second reduced GGQ rule.");

    /// Export enum GridType
    py::enum_<GridType>(m, "GridType")
        .value("B_SPLINE_GRID", GridType::b_spline_grid, "B-spline background grid.")
        .value("HEXAHEDRAL_FE_GRID", GridType::hexahedral_fe_grid, "Hexahedral finite-element grid.");
}

}  // namespace queso::python

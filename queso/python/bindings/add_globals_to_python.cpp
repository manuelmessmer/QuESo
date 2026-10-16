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

// External includes
#include <pybind11/native_enum.h>

// Project includes
#include "queso/python/bindings/add_globals_to_python.h"
// To export
#include "queso/includes/define.hpp"
#include "queso/utilities/triangle_gauss_legendre_integration_points.hpp"

namespace queso::python {

namespace py = pybind11;

void AddGlobalsToPython(pybind11::module& m)
{
    using triangle_utilities::TriangleQuadratureOrder;
    py::native_enum<TriangleQuadratureOrder>(
        m,
        "TriangleQuadratureOrder",
        "enum.Enum",
        "Maximum polynomial degree integrated exactly by triangle quadrature.\n\n"
        "Members:\n"
        "- DEGREE_1: One point. Exact through degree 1.\n"
        "- DEGREE_2: Three points. Exact through degree 2.\n"
        "- DEGREE_4: Six points. Exact through degree 4.\n"
        "- DEGREE_6: Twelve points. Exact through degree 6."
    )
        .value("DEGREE_1", TriangleQuadratureOrder::Degree1)
        .value("DEGREE_2", TriangleQuadratureOrder::Degree2)
        .value("DEGREE_4", TriangleQuadratureOrder::Degree4)
        .value("DEGREE_6", TriangleQuadratureOrder::Degree6)
        .finalize();

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

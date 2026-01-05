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
#include "queso/python/bindings/define_python.hpp"
#include "queso/python/bindings/add_globals_to_python.h"
// To export
#include "queso/includes/define.hpp"

namespace queso {
namespace Python {

namespace py = pybind11;

void AddGlobalsToPython(pybind11::module& m) {

    /// Export enum IntegrationMethod
    py::enum_<IntegrationMethod>(m, "IntegrationMethod")
        .value("Gauss", IntegrationMethod::gauss)
        .value("Gauss_Reduced1", IntegrationMethod::gauss_reduced_1)
        .value("Gauss_Reduced2", IntegrationMethod::gauss_reduced_2)
        .value("GGQ_Optimal", IntegrationMethod::ggq_optimal)
        .value("GGQ_Reduced1", IntegrationMethod::ggq_reduced_1)
        .value("GGQ_Reduced2", IntegrationMethod::ggq_reduced_2)
    ;

    /// Export enum GridType
    py::enum_<GridType>(m, "GridType")
        .value("b_spline_grid", GridType::b_spline_grid)
        .value("hexahedral_fe_grid", GridType::hexahedral_fe_grid)
    ;

} // End AddGlobalsToPython

} // End namespace Python
} // End namespace queso


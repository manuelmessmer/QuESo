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
#include "queso/python/define_python.hpp"
#include "queso/python/add_globals_to_python.h"
// To export
#include "queso/includes/define.hpp"

namespace queso {
namespace Python {

namespace py = pybind11;

void AddGlobalsToPython(pybind11::module& m) {

    /// Export enum IntegrationMethod
    py::enum_<IntegrationMethod>(m, "IntegrationMethod")
        .value("Gauss", IntegrationMethod::Gauss)
        .value("Gauss_Reduced1", IntegrationMethod::Gauss_Reduced1)
        .value("Gauss_Reduced2", IntegrationMethod::Gauss_Reduced2)
        .value("GGQ_Optimal", IntegrationMethod::GGQ_Optimal)
        .value("GGQ_Reduced1", IntegrationMethod::GGQ_Reduced1)
        .value("GGQ_Reduced2", IntegrationMethod::GGQ_Reduced2)
        .export_values()
    ;

    /// Export enum GridType
    py::enum_<GridType>(m, "GridType")
        .value("b_spline_grid", GridType::b_spline_grid)
        .value("hexahedral_fe_grid", GridType::hexahedral_fe_grid)
    ;

} // End AddGlobalsToPython

} // End namespace Python
} // End namespace queso


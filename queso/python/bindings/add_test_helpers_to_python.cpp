//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

// External includes
#include <pybind11/stl.h>

// Project includes
#include "queso/python/bindings/add_test_helpers_to_python.h"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso::python {

namespace py = pybind11;

void AddTestHelpersToPython(py::module& rModule)
{
    using namespace py::literals;

    py::class_<IntegrationPointFactory1D>(rModule, "IntegrationPointFactory1D")
        .def_static(
            "get_ggq",
            [](SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethodType Method) {
                auto pPoints = IntegrationPointFactory1D::GetGGQ(PolynomialDegree, NumberKnotSpans, Method);
                return std::move(*pPoints);
            },
            "polynomial_degree"_a,
            "num_elements"_a,
            "integration_method"_a,
            "Return a generalized Gaussian quadrature rule as ``[(coordinate, weight), ...]``."
        );
}

}  // namespace queso::python

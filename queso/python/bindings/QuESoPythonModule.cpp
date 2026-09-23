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
#include <pybind11/pybind11.h>

// Project includes
#include "queso/includes/define.hpp"

// To export
#include "queso/python/bindings/add_containers_to_python.h"
#include "queso/python/bindings/add_dictionary_to_python.h"
#include "queso/python/bindings/add_globals_to_python.h"
#include "queso/python/bindings/add_io_to_python.h"
#if defined(QUESO_PYTHON_BUILD_TEST_HELPERS)
#include "queso/python/bindings/add_test_helpers_to_python.h"
#endif
#include "queso/python/bindings/add_utilities_to_python.h"

namespace queso::python {

namespace py = pybind11;

namespace {
    void AddAllToPython(py::module& rPyQuESo)
    {
        auto MeshModule = rPyQuESo.def_submodule("mesh", "Triangle mesh types and algorithms.");
        auto IoModule = rPyQuESo.def_submodule("io", "Mesh and settings input/output functions.");

        AddGlobalsToPython(rPyQuESo);
        AddDictionaryToPython(rPyQuESo);
        AddContainersToPython(rPyQuESo, MeshModule);
        AddUtilitiesToPython(MeshModule);
        AddIoToPython(IoModule);

#if defined(QUESO_PYTHON_BUILD_TEST_HELPERS)
        auto TestingModule = rPyQuESo.def_submodule("testing", "Testing-only helpers.");
        AddTestHelpersToPython(TestingModule);
#endif
    }
}  // namespace

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Private extension module for pyqueso.";

    m.def(
        "print_logo",
        []() {
            QuESo_INFO << " Importing QuESo \n"
                       << "   ____        ______  _____        \n"
                       << "  / __ \\      |  ____|/ ____|       \n"
                       << " | |  | |_   _| |__  | (___   ___   \n"
                       << " | |  | | | | |  __|  \\___ \\ / _ \\  \n"
                       << " | |__| | |_| | |____ ____) | (_) | \n"
                       << "  \\___\\_\\\\__,_|______|_____/ \\___/  \n"
                       << "\t Quadrature for Embedded Solids \n\n";
        },
        "Print the QuESo logo."
    );

    auto Sys = py::module::import("sys");
    auto PyQuESo = Sys.attr("modules")["pyqueso"].cast<py::module>();
    AddAllToPython(PyQuESo);
}

}  // namespace queso::python

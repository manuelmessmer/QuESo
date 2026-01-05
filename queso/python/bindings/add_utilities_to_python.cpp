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
#include "queso/python/bindings/add_utilities_to_python.h"

/// To export
#include "queso/utilities/mesh_utilities.h"

namespace queso {
namespace Python {

namespace py = pybind11;

void AddUtilitiesToPython(pybind11::module& m) {

    /// Export MeshUtilites
    py::class_<MeshUtilities>(m,"MeshUtilities")
        .def_static("Volume", &MeshUtilities::VolumeOMP)
        .def_static("Area", &MeshUtilities::AreaOMP)
    ;

} // End AddUtilitiesToPython

} // End namespace Python
} // End namespace queso


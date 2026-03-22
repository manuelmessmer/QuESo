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
    auto mesh_util = m.def_submodule("MeshUtil");

    mesh_util.def("Volume", [](const TriangleMesh &rMesh) {
        return MeshUtilities::VolumeOMP(rMesh.View());
    });
    mesh_util.def("Area", [](const TriangleMesh &rMesh) {
        return MeshUtilities::AreaOMP(rMesh.View());
    });

} // End AddUtilitiesToPython

} // End namespace Python
} // End namespace queso

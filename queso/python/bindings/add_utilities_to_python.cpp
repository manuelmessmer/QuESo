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

namespace queso::python {

namespace py = pybind11;

void AddUtilitiesToPython(pybind11::module& m)
{
    m.def(
        "volume",
        [](const TriangleMesh& rMesh) { return MeshUtilities::VolumeOMP(rMesh.View()); },
        py::arg("mesh"),
        "Compute the signed volume enclosed by a triangle mesh."
    );
    m.def(
        "area",
        [](const TriangleMesh& rMesh) { return MeshUtilities::AreaOMP(rMesh.View()); },
        py::arg("mesh"),
        "Compute the surface area of a triangle mesh."
    );
}

}  // namespace queso::python

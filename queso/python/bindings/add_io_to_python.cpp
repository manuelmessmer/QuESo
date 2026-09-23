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
#include "queso/python/bindings/add_io_to_python.h"

/// To export
#include "queso/io/io_utilities.h"

namespace queso::python {

namespace py = pybind11;

void AddIoToPython(pybind11::module& m)
{
    using DictionaryType = Dictionary<key::MainValuesTypeTag>;

    m.def(
        "read_mesh_from_stl",
        &IO::ReadMeshFromSTL,
        py::arg("mesh"),
        py::arg("filename"),
        "Read an STL file into a mutable triangle mesh."
    );
    m.def(
        "write_dictionary_to_json",
        &IO::WriteDictionaryToJSON<DictionaryType>,
        py::arg("dictionary"),
        py::arg("filename"),
        "Write a dictionary to a JSON file."
    );
}

}  // namespace queso::python

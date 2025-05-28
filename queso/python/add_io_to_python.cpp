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
#include "queso/python/add_io_to_python.h"
// To export
#include "queso/io/io_utilities.h"

namespace queso {
namespace Python {

namespace py = pybind11;

void AddIoToPython(pybind11::module& m) {
    using DictionaryType = Dictionary<key::MainValuesTypeTag>;

    /// Export IoUtilities
    py::class_<IO>(m,"IO")
        .def_static("ReadMeshFromSTL", &IO::ReadMeshFromSTL)
        .def_static("WriteDictionaryToJSON", &IO::WriteDictionaryToJSON<DictionaryType>)
    ;

} // End AddIoToPython

} // End namespace Python
} // End namespace queso


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

#ifndef ADD_IO_TO_PYTHON_INCLUDE_H
#define ADD_IO_TO_PYTHON_INCLUDE_H

// External includes
#include <pybind11/pybind11.h>

namespace queso {
namespace Python {

    void AddIoToPython(pybind11::module& m);

} // End namespace Python
} // End namespace queso

#endif // ADD_IO_TO_PYTHON_INCLUDE_H
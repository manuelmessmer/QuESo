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

#ifndef DEFINE_PYTHON_INCLUDE_HPP
#define DEFINE_PYTHON_INCLUDE_HPP

/// External includes
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace queso {
namespace Python {

template<class TType>
std::string PrintObject(const TType& rObject)
{
    std::stringstream ss;
    ss << rObject;
    return ss.str();
}


} // End namespace Python
} // End namespace queso

#endif // DEFINE_PYTHON_INCLUDE_HPP
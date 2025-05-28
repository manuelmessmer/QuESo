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
#include <pybind11/functional.h>

/// Project includes
#include "queso/includes/define.hpp"

namespace queso {
namespace Python {

/// @brief Allows to pass unique ownership from Python to C++.
/// @tparam TType
template<typename TType>
struct UniqueHolder {
    /// @brief Constructor (Takes ownership).
    /// @param NewData
    UniqueHolder(Unique<TType> NewData) : mpData(std::move(NewData))
    {
    }

    /// @brief Returns reference to the underlying object.
    /// @return TType&
    TType& GetObject(){
        QuESo_ERROR_IF(!mpData) << "The underlying data has been released already.\n";
        return *mpData;
    }

    /// @brief Releases the ownership of the underlying object.
    /// @return Unique<TType>
    Unique<TType> Release() {
        QuESo_ERROR_IF(!mpData) << "The underlying data has been released already.\n";
        return std::move(mpData);
    }

    /// Ptr to object.
    Unique<TType> mpData;
};

/// @brief Helper function to print an object to a string using stringstream.
/// @tparam TType
/// @param rObject
/// @return std::string
template<class TType>
std::string PrintObject(const TType& rObject) {
    std::stringstream ss;
    ss << rObject;
    return ss.str();
}

} // End namespace Python
} // End namespace queso

#endif // DEFINE_PYTHON_INCLUDE_HPP
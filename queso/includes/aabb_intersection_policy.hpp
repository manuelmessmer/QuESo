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

#pragma once

//// STL includes
#include <cstdint>

namespace queso {

/// @brief Selects whether AABB boundary touching counts as an intersection.
enum class AabbIntersectionPolicy : std::uint8_t {
    Exact,  ///< Tests the original AABB, so touching geometry intersects.
    SnapEroded  ///< Erodes the AABB by the resolved snap distance before testing.
};

}  // namespace queso

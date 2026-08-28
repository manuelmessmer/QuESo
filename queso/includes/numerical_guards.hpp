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
#include <cmath>

namespace queso::numerical_guards {

/// Default absolute tolerance used by generic numerical guards.
inline constexpr double Small = 1e-14;

/// @brief Returns whether a finite value is safely above a minimum divisor magnitude.
/// @param Value Candidate divisor.
/// @param MinimumMagnitude Finite non-negative magnitude required for a safe division.
/// @return True when Value is finite and its magnitude exceeds MinimumMagnitude.
[[nodiscard]] inline bool IsSafeDivisor(double Value, double MinimumMagnitude = Small) noexcept
{
    return std::isfinite(Value) && std::isfinite(MinimumMagnitude) && MinimumMagnitude >= 0.0
           && std::abs(Value) > MinimumMagnitude;
}

}  // namespace queso::numerical_guards

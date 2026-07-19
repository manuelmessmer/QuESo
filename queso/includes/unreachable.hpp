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

#include "queso/includes/assert.hpp"

namespace queso {

/// @brief Indicates that an unreachable code path has been executed.
/// @param Msg Message displayed by the assertion in debug builds.
/// @note This function is intended only for code paths that are unreachable
///       in a correct program. Reaching it results in undefined behavior in
///       release builds.
[[noreturn]] inline void
    Unreachable([[maybe_unused]] std::string_view Msg = "Reached an unreachable code path") noexcept(NOTDEBUG)
{
    QuESo_ASSERT(false, Msg);

#if defined(__GNUC__) || defined(__clang__)
    __builtin_unreachable();
#elif defined(_MSC_VER)
    __assume(false);
#else
    std::abort();
#endif
}

}  // namespace queso

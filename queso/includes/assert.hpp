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

/// STL includes
#include <assert.h>
#include <type_traits>

/// Project includes
// IWYU pragma: begin_exports
#include "queso/includes/exception.hpp"
// IWYU pragma: end_exports

namespace queso {

/// QuESo Macros
#ifdef NDEBUG  // asserts disabled

constexpr bool NOTDEBUG = true;

#define QuESo_ASSERT(Assertation, Message) (assert(0))

#else  // asserts enabled

constexpr bool NOTDEBUG = false;

namespace detail {
    template<typename A>
    constexpr void
        Assert(A assertion, std::string_view message, const char* file, const char* function, std::size_t line)
    {
        if (assertion) return;

        if (std::is_constant_evaluated()) {
            throw message;
        } else {
            throw Exception(file, function, line) << message;
        }
    }

#define QuESo_ASSERT(Assertation, Message) \
    queso::detail::Assert(Assertation, Message, __FILE__, QuESo_CURRENT_FUNCTION, __LINE__)

}  // namespace detail

#endif  // NDEBUG

}  // namespace queso

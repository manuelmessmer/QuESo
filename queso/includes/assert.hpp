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
#include <source_location>
#include <string_view>

/// Project includes
// IWYU pragma: begin_exports
#include "queso/includes/exception.hpp"
// IWYU pragma: end_exports

namespace queso {

/// QuESo Macros
#ifdef NDEBUG  // asserts disabled

constexpr bool NOTDEBUG = true;

#else  // asserts enabled

constexpr bool NOTDEBUG = false;

#endif  // NDEBUG

namespace detail {
    /// @brief Reports an assertion failure at the caller's source location.
    /// @param Expression Failed expression text.
    /// @param Message Explanation of the violated requirement.
    /// @param Location Assertion call site.
    [[noreturn]] inline void AssertionFailed(
        std::string_view Expression,
        std::string_view Message,
        std::source_location Location = std::source_location::current()
    )
    {
        throw Exception(Location.file_name(), Location.function_name(), Location.line())
            << "Assertion failed: " << Expression << ". " << Message;
    }

    /// @brief Returns whether optional validation should execute.
    /// @details The release environment setting is cached on first use in the library.
    /// @return Always true in debug builds. In release, true only when QUESO_DEBUG was exactly "1" at initialization.
    [[nodiscard]] bool DiagnosticsEnabled() noexcept;
}  // namespace detail

// Conditions are evaluated once. Messages are evaluated only on failure.
// Disabled assertions evaluate neither argument. Conditions must have no side effects.
#define QuESo_CONTRACT_ASSERT(Condition, Message)                                      \
    do {                                                                               \
        if (!(Condition)) { ::queso::detail::AssertionFailed(#Condition, (Message)); } \
    } while (false)

#ifndef NDEBUG
#define QuESo_DEBUG_ASSERT(Condition, Message) QuESo_CONTRACT_ASSERT(Condition, Message)
#else
#define QuESo_DEBUG_ASSERT(Condition, Message) ((void)0)
#endif

#define QuESo_DIAGNOSTICS_ENABLED() (::queso::detail::DiagnosticsEnabled())
#define QuESo_DIAGNOSTIC_ASSERT(Condition, Message)                                     \
    do {                                                                                \
        if (QuESo_DIAGNOSTICS_ENABLED()) { QuESo_CONTRACT_ASSERT(Condition, Message); } \
    } while (false)

}  // namespace queso

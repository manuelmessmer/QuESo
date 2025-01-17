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

#ifndef ASSERT_HPP_INCLUDE
#define ASSERT_HPP_INCLUDE

/// STL includes
#include <assert.h>
/// Project includes
#include "queso/includes/define.hpp"
#include "queso/includes/exception.hpp"

namespace queso {

/// QuESo Macros
#ifdef NDEBUG    // asserts disabled

static constexpr bool NOTDEBUG = true;
// Use empty assert from STL (Different for Win and GCC)
# define QuESo_ASSERT(Assertation, Message)		(assert(0))

#else            // asserts enabled

static constexpr bool NOTDEBUG = false;

namespace detail {

template<typename A>
inline void Assert(A assertion, const std::string& rMessage,
                   std::string const& rFileName, std::string const& rFunctionName, std::size_t LineNumber) {
    if( !assertion ) {
        throw Exception(rFileName, rFunctionName, LineNumber) << rMessage;
    }
} // End namespace detail

#define QuESo_ASSERT(Assertation, Message) queso::detail::Assert(Assertation, Message, __FILE__, QuESo_CURRENT_FUNCTION, __LINE__)

} // End detail

#endif // NDEBUG

} // End queso

#endif // ASSERT_HPP_INCLUDE
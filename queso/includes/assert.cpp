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

//// STL includes
#include <cstdlib>
#include <string_view>

//// Project includes
#include "queso/includes/assert.hpp"

namespace queso::detail {

bool DiagnosticsEnabled() noexcept
{
#ifndef NDEBUG
    return true;
#else
    static const bool is_enabled = []() noexcept {
        const char* p_value = std::getenv("QUESO_DEBUG");
        return p_value != nullptr && std::string_view(p_value) == "1";
    }();
    return is_enabled;
#endif
}

}  // namespace queso::detail

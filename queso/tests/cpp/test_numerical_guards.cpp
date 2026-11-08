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

//// External includes
#include <doctest/doctest.h>

//// STL includes
#include <limits>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/includes/numerical_guards.hpp"

namespace queso::testing {

TEST_SUITE("NumericalGuardsTestSuite")
{

    TEST_CASE("NumericalGuardsTestSuite - Safe Divisor Requires Finite Magnitude Above Threshold")
    {
        CHECK_FALSE(numerical_guards::IsSafeDivisor(0.0));
        CHECK_FALSE(numerical_guards::IsSafeDivisor(numerical_guards::Small));
        CHECK(numerical_guards::IsSafeDivisor(2.0 * numerical_guards::Small));
        CHECK_FALSE(numerical_guards::IsSafeDivisor(std::numeric_limits<double>::infinity()));
        CHECK_FALSE(numerical_guards::IsSafeDivisor(std::numeric_limits<double>::quiet_NaN()));
    }

}  // TEST_SUITE

}  // namespace queso::testing

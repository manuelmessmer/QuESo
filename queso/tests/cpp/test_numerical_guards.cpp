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
#include <boost/test/unit_test.hpp>

//// STL includes
#include <limits>

//// Project includes
#include "queso/includes/checks.hpp"
#include "queso/includes/numerical_guards.hpp"

namespace queso::Testing {

BOOST_AUTO_TEST_SUITE(NumericalGuardsTestSuite)

BOOST_AUTO_TEST_CASE(SafeDivisorRequiresFiniteMagnitudeAboveThreshold)
{
    QuESo_CHECK_IS_FALSE(numerical_guards::IsSafeDivisor(0.0));
    QuESo_CHECK_IS_FALSE(numerical_guards::IsSafeDivisor(numerical_guards::Small));
    QuESo_CHECK(numerical_guards::IsSafeDivisor(2.0 * numerical_guards::Small));
    QuESo_CHECK_IS_FALSE(numerical_guards::IsSafeDivisor(std::numeric_limits<double>::infinity()));
    QuESo_CHECK_IS_FALSE(numerical_guards::IsSafeDivisor(std::numeric_limits<double>::quiet_NaN()));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

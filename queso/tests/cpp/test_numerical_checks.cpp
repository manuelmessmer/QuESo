// Numerical comparison contracts used by the geometry tests.
#include "helpers/numerical_checks.hpp"
#include <array>
#include <doctest/doctest.h>
#include <limits>

namespace queso::testing {
TEST_SUITE("NumericalChecks")
{

    TEST_CASE("NumericalChecks - absolute tolerance is inclusive")
    {
        CHECK(IsNear(1.0, 1.25, 0.25));
        CHECK_FALSE(IsNear(1.0, 1.25, 0.125));
        CHECK(IsNear(-1.0, -1.25, 0.25));
        CHECK(IsNear(1.0, 1.0, 0.0));
    }

    TEST_CASE("NumericalChecks - relative tolerance uses the expected value")
    {
        CHECK(IsRelativeNear(1.5, 1.0, 0.5));
        CHECK_FALSE(IsRelativeNear(1.5, 1.0, 0.4));
        CHECK(IsRelativeNear(-1.5, -1.0, 0.5));
        CHECK(IsRelativeNear(0.0009, 0.0, 0.001));
        CHECK_FALSE(IsRelativeNear(0.0009, 0.0, 0.0001));
        const double epsilon = std::numeric_limits<double>::epsilon();
        CHECK(IsRelativeNear(epsilon + 0.25, epsilon, 0.25));
    }

    TEST_CASE("NumericalChecks - absolute point comparisons accept all coordinates within tolerance")
    {
        constexpr std::array<double, 3> expected{ 2.0, -4.0, 0.0 };
        constexpr std::array<double, 3> actual{ 2.125, -4.125, 0.125 };
        CheckPointNear(actual, expected, 0.125);
        CheckPointNear(expected, expected, 0.0);
    }

    TEST_CASE(
        "NumericalChecks - absolute point comparisons reject a mismatch in each coordinate"
        * doctest::expected_failures(3) * doctest::no_breaks()
    )
    {
        constexpr std::array<double, 3> expected{ 2.0, -4.0, 0.0 };
        // Exactly one helper assertion must fail for each independently perturbed coordinate.
        for (unsigned axis = 0; axis < expected.size(); ++axis) {
            CAPTURE(axis);
            auto actual = expected;
            actual[axis] += 0.25;
            CheckPointNear(actual, expected, 0.125);
        }
    }

    TEST_CASE("NumericalChecks - relative point comparisons preserve the near-zero fallback")
    {
        constexpr std::array<double, 3> expected{ 2.0, -4.0, 0.0 };
        constexpr std::array<double, 3> actual{ 2.25, -4.5, 0.125 };
        CheckPointRelativeNear(actual, expected, 0.125);
        CheckPointRelativeNear(expected, expected, 0.0);
    }

    TEST_CASE(
        "NumericalChecks - relative point comparisons reject a mismatch in each coordinate"
        * doctest::expected_failures(3) * doctest::no_breaks()
    )
    {
        constexpr std::array<double, 3> expected{ 2.0, -4.0, 0.0 };
        constexpr std::array<double, 3> offsets{ 0.5, -1.0, 0.25 };
        // Covers positive/negative reference values and the absolute fallback at zero.
        for (unsigned axis = 0; axis < expected.size(); ++axis) {
            CAPTURE(axis);
            auto actual = expected;
            actual[axis] += offsets[axis];
            CheckPointRelativeNear(actual, expected, 0.125);
        }
    }

    TEST_CASE("NumericalChecks - NaN never satisfies a tolerance")
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(IsNear(nan, 1.0, 1.0));
        CHECK_FALSE(IsRelativeNear(1.0, nan, 1.0));
    }

}  // TEST_SUITE
}  // namespace queso::testing

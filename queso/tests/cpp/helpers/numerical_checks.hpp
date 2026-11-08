#pragma once

//// STL includes
#include <cmath>
#include <limits>
#include <source_location>

//// External includes
#include <doctest/doctest.h>

namespace queso::testing {

/// Returns whether two values satisfy an inclusive absolute tolerance.
[[nodiscard]] inline bool IsNear(double Actual, double Expected, double Tolerance) noexcept
{ return std::abs(Actual - Expected) <= Tolerance; }

/// Preserves the reference-relative tolerance and absolute near-zero fallback.
[[nodiscard]] inline bool IsRelativeNear(double Actual, double Expected, double Tolerance) noexcept
{
    return std::abs(Expected) <= std::numeric_limits<double>::epsilon()
               ? IsNear(Actual, Expected, Tolerance)
               : std::abs((Actual - Expected) / Expected) <= Tolerance;
}

/// Checks a scalar using inclusive absolute tolerance with failure and caller context.
inline void CheckNear(
    double Actual,
    double Expected,
    double Tolerance,
    std::source_location Location = std::source_location::current()
)
{
    INFO("comparison requested at ", doctest::String(Location.file_name()), ":", Location.line());
    INFO("actual: ", Actual, ", expected: ", Expected, ", absolute tolerance: ", Tolerance);
    CHECK(IsNear(Actual, Expected, Tolerance));
}

/// Checks a scalar using reference-relative tolerance and an absolute near-zero fallback.
inline void CheckRelativeNear(
    double Actual,
    double Expected,
    double Tolerance,
    std::source_location Location = std::source_location::current()
)
{
    INFO("comparison requested at ", doctest::String(Location.file_name()), ":", Location.line());
    INFO("actual: ", Actual, ", expected: ", Expected, ", relative tolerance: ", Tolerance);
    CHECK(IsRelativeNear(Actual, Expected, Tolerance));
}

/// Checks each coordinate independently with failure context.
template<typename TActual, typename TExpected>
void CheckPointNear(
    const TActual& rActual,
    const TExpected& rExpected,
    double Tolerance,
    std::source_location Location = std::source_location::current()
)
{
    INFO("comparison requested at ", doctest::String(Location.file_name()), ":", Location.line());
    for (unsigned axis = 0; axis < 3; ++axis) {
        CAPTURE(axis);
        INFO("actual: ", rActual[axis], ", expected: ", rExpected[axis], ", absolute tolerance: ", Tolerance);
        CHECK(IsNear(rActual[axis], rExpected[axis], Tolerance));
    }
}

/// Checks each coordinate using the existing reference-relative tolerance.
template<typename TActual, typename TExpected>
void CheckPointRelativeNear(
    const TActual& rActual,
    const TExpected& rExpected,
    double Tolerance,
    std::source_location Location = std::source_location::current()
)
{
    INFO("comparison requested at ", doctest::String(Location.file_name()), ":", Location.line());
    for (unsigned axis = 0; axis < 3; ++axis) {
        CAPTURE(axis);
        INFO("actual: ", rActual[axis], ", expected: ", rExpected[axis], ", relative tolerance: ", Tolerance);
        CHECK(IsRelativeNear(rActual[axis], rExpected[axis], Tolerance));
    }
}

}  // namespace queso::testing

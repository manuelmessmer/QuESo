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
#include "queso/containers/geometry_tolerance.hpp"

namespace queso::testing {

TEST_SUITE("GeometryToleranceTestSuite")
{

    TEST_CASE("GeometryToleranceTestSuite - Resolves Dimensionally At Unit Scale")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

        CHECK((tolerance.SnapDistance()) > (tolerance.ZeroLength()));
        CHECK((tolerance.ZeroArea()) == (tolerance.ZeroLength() * tolerance.ZeroLength()));
        CHECK((tolerance.ZeroVolume()) == (tolerance.ZeroArea() * tolerance.ZeroLength()));
    }

    TEST_CASE("GeometryToleranceTestSuite - Uses Coordinate Resolution At Small And Translated Scales")
    {
        constexpr double roundoff_multiplier = 32.0;
        constexpr double small_length = 1e-6;
        const double unit_resolution = roundoff_multiplier * std::numeric_limits<double>::epsilon();
        const GeometryTolerance small =
            GeometryTolerance::FromScale({ .length_scale = small_length, .coordinate_scale = 1.0 });

        CHECK((small.SnapDistance()) == (unit_resolution));
        CHECK((small.ZeroLength()) == (unit_resolution));
        CHECK((small.ZeroArea()) == (unit_resolution * unit_resolution));
        CHECK((small.ZeroVolume()) == (unit_resolution * unit_resolution * unit_resolution));

        constexpr double coordinate_scale = 1e12;
        const double translated_resolution = unit_resolution * coordinate_scale;
        const GeometryTolerance translated =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = coordinate_scale });

        CHECK((translated.SnapDistance()) == (translated_resolution));
        CHECK((translated.ZeroLength()) == (translated_resolution));
        CHECK((translated.ZeroArea()) == (translated_resolution * translated_resolution));
        CHECK((translated.ZeroVolume()) == (translated_resolution * translated_resolution * translated_resolution));
    }

    TEST_CASE("GeometryToleranceTestSuite - Uses Relative Thresholds At Large Scale")
    {
        constexpr double length = 1e6;
        const GeometryTolerance unit = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = length, .coordinate_scale = length });

        CheckRelativeNear(tolerance.SnapDistance(), unit.SnapDistance() * length, EPS4);
        CheckRelativeNear(tolerance.ZeroLength(), unit.ZeroLength() * length, EPS4);
        CheckRelativeNear(tolerance.ZeroArea(), unit.ZeroArea() * length * length, EPS4);
        CheckRelativeNear(tolerance.ZeroVolume(), unit.ZeroVolume() * length * length * length, EPS4);
    }

    TEST_CASE("GeometryToleranceTestSuite - Resolves Scales From Bounds")
    {
        const BoundingBoxType bounds = MakeBox({ -3.0, -2.0, -1.0 }, { 5.0, 1.0, 2.0 });
        const GeometryTolerance from_bounds = GeometryTolerance::FromBounds(bounds);
        const GeometryTolerance from_scale =
            GeometryTolerance::FromScale({ .length_scale = 8.0, .coordinate_scale = 5.0 });

        CHECK((from_bounds.SnapDistance()) == (from_scale.SnapDistance()));
        CHECK((from_bounds.ZeroLength()) == (from_scale.ZeroLength()));
        CHECK((from_bounds.ZeroArea()) == (from_scale.ZeroArea()));
        CHECK((from_bounds.ZeroVolume()) == (from_scale.ZeroVolume()));
    }

    TEST_CASE("GeometryToleranceTestSuite - Applies Inclusive Dimensional Predicates")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const PointType origin{};
        const PointType same_point{ tolerance.SnapDistance(), 0.0, 0.0 };
        const PointType different_point{ 2.0 * tolerance.SnapDistance(), 0.0, 0.0 };

        CHECK(tolerance.CoordinatesAreSame(0.0, tolerance.SnapDistance()));
        CHECK_FALSE(tolerance.CoordinatesAreSame(0.0, 2.0 * tolerance.SnapDistance()));
        CHECK(tolerance.PointsAreSame(origin, same_point));
        CHECK_FALSE(tolerance.PointsAreSame(origin, different_point));
        CHECK(tolerance.IsZeroLength(-tolerance.ZeroLength()));
        CHECK_FALSE(tolerance.IsZeroLength(2.0 * tolerance.ZeroLength()));
        CHECK(tolerance.IsZeroArea(-tolerance.ZeroArea()));
        CHECK_FALSE(tolerance.IsZeroArea(2.0 * tolerance.ZeroArea()));
        CHECK(tolerance.IsZeroVolume(-tolerance.ZeroVolume()));
        CHECK_FALSE(tolerance.IsZeroVolume(2.0 * tolerance.ZeroVolume()));
    }

    TEST_CASE("GeometryToleranceTestSuite - Rejects Invalid Or Unrepresentable Scales")
    {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        constexpr double infinity = std::numeric_limits<double>::infinity();

        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 0.0, .coordinate_scale = 1.0 })), Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = -1.0, .coordinate_scale = 1.0 })),
            Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = nan, .coordinate_scale = 1.0 })), Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 0.0 })), Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = infinity })),
            Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1e-16, .coordinate_scale = 1.0 })),
            Exception
        );
        CHECK_THROWS_AS(
            static_cast<void>(GeometryTolerance::FromBounds(MakeBox({ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 1.0 }))), Exception
        );
    }

}  // TEST_SUITE

}  // namespace queso::testing

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
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/includes/checks.hpp"

namespace queso::Testing {

BOOST_AUTO_TEST_SUITE(GeometryToleranceTestSuite)

BOOST_AUTO_TEST_CASE(ResolvesDimensionallyAtUnitScale)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

    QuESo_CHECK_GT(tolerance.SnapDistance(), tolerance.ZeroLength());
    QuESo_CHECK_EQUAL(tolerance.ZeroArea(), tolerance.ZeroLength() * tolerance.ZeroLength());
    QuESo_CHECK_EQUAL(tolerance.ZeroVolume(), tolerance.ZeroArea() * tolerance.ZeroLength());
}

BOOST_AUTO_TEST_CASE(UsesCoordinateResolutionAtSmallAndTranslatedScales)
{
    constexpr double roundoff_multiplier = 32.0;
    constexpr double small_length = 1e-6;
    const double unit_resolution = roundoff_multiplier * std::numeric_limits<double>::epsilon();
    const GeometryTolerance small =
        GeometryTolerance::FromScale({ .length_scale = small_length, .coordinate_scale = 1.0 });

    QuESo_CHECK_EQUAL(small.SnapDistance(), unit_resolution);
    QuESo_CHECK_EQUAL(small.ZeroLength(), unit_resolution);
    QuESo_CHECK_EQUAL(small.ZeroArea(), unit_resolution * unit_resolution);
    QuESo_CHECK_EQUAL(small.ZeroVolume(), unit_resolution * unit_resolution * unit_resolution);

    constexpr double coordinate_scale = 1e12;
    const double translated_resolution = unit_resolution * coordinate_scale;
    const GeometryTolerance translated =
        GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = coordinate_scale });

    QuESo_CHECK_EQUAL(translated.SnapDistance(), translated_resolution);
    QuESo_CHECK_EQUAL(translated.ZeroLength(), translated_resolution);
    QuESo_CHECK_EQUAL(translated.ZeroArea(), translated_resolution * translated_resolution);
    QuESo_CHECK_EQUAL(translated.ZeroVolume(), translated_resolution * translated_resolution * translated_resolution);
}

BOOST_AUTO_TEST_CASE(UsesRelativeThresholdsAtLargeScale)
{
    constexpr double length = 1e6;
    const GeometryTolerance unit = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const GeometryTolerance tolerance =
        GeometryTolerance::FromScale({ .length_scale = length, .coordinate_scale = length });

    QuESo_CHECK_RELATIVE_NEAR(tolerance.SnapDistance(), unit.SnapDistance() * length, EPS4);
    QuESo_CHECK_RELATIVE_NEAR(tolerance.ZeroLength(), unit.ZeroLength() * length, EPS4);
    QuESo_CHECK_RELATIVE_NEAR(tolerance.ZeroArea(), unit.ZeroArea() * length * length, EPS4);
    QuESo_CHECK_RELATIVE_NEAR(tolerance.ZeroVolume(), unit.ZeroVolume() * length * length * length, EPS4);
}

BOOST_AUTO_TEST_CASE(ResolvesScalesFromBounds)
{
    const BoundingBoxType bounds = MakeBox({ -3.0, -2.0, -1.0 }, { 5.0, 1.0, 2.0 });
    const GeometryTolerance from_bounds = GeometryTolerance::FromBounds(bounds);
    const GeometryTolerance from_scale = GeometryTolerance::FromScale({ .length_scale = 8.0, .coordinate_scale = 5.0 });

    QuESo_CHECK_EQUAL(from_bounds.SnapDistance(), from_scale.SnapDistance());
    QuESo_CHECK_EQUAL(from_bounds.ZeroLength(), from_scale.ZeroLength());
    QuESo_CHECK_EQUAL(from_bounds.ZeroArea(), from_scale.ZeroArea());
    QuESo_CHECK_EQUAL(from_bounds.ZeroVolume(), from_scale.ZeroVolume());
}

BOOST_AUTO_TEST_CASE(AppliesInclusiveDimensionalPredicates)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const PointType origin{};
    const PointType same_point{ tolerance.SnapDistance(), 0.0, 0.0 };
    const PointType different_point{ 2.0 * tolerance.SnapDistance(), 0.0, 0.0 };

    QuESo_CHECK(tolerance.CoordinatesAreSame(0.0, tolerance.SnapDistance()));
    QuESo_CHECK_IS_FALSE(tolerance.CoordinatesAreSame(0.0, 2.0 * tolerance.SnapDistance()));
    QuESo_CHECK(tolerance.PointsAreSame(origin, same_point));
    QuESo_CHECK_IS_FALSE(tolerance.PointsAreSame(origin, different_point));
    QuESo_CHECK(tolerance.IsZeroLength(-tolerance.ZeroLength()));
    QuESo_CHECK_IS_FALSE(tolerance.IsZeroLength(2.0 * tolerance.ZeroLength()));
    QuESo_CHECK(tolerance.IsZeroArea(-tolerance.ZeroArea()));
    QuESo_CHECK_IS_FALSE(tolerance.IsZeroArea(2.0 * tolerance.ZeroArea()));
    QuESo_CHECK(tolerance.IsZeroVolume(-tolerance.ZeroVolume()));
    QuESo_CHECK_IS_FALSE(tolerance.IsZeroVolume(2.0 * tolerance.ZeroVolume()));
}

BOOST_AUTO_TEST_CASE(RejectsInvalidOrUnrepresentableScales)
{
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double infinity = std::numeric_limits<double>::infinity();

    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 0.0, .coordinate_scale = 1.0 })), Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = -1.0, .coordinate_scale = 1.0 })), Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = nan, .coordinate_scale = 1.0 })), Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 0.0 })), Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = infinity })),
        Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromScale({ .length_scale = 1e-16, .coordinate_scale = 1.0 })), Exception
    );
    BOOST_CHECK_THROW(
        static_cast<void>(GeometryTolerance::FromBounds(MakeBox({ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 1.0 }))), Exception
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

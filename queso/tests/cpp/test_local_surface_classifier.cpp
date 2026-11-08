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
#include <array>

//// External includes
#include <doctest/doctest.h>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/embedding/local_surface_classifier.h"

namespace queso::testing {
namespace {

    [[nodiscard]] TriangleMesh MakeOrientedSection(double Scale)
    {
        TriangleMesh mesh;
        const IndexType first = mesh.AddVertex({ Scale, 0.0, 0.0 });
        const IndexType second = mesh.AddVertex({ Scale, Scale, 0.0 });
        const IndexType third = mesh.AddVertex({ Scale, 0.0, Scale });
        mesh.AddTriangle({ first, second, third }, { 1.0, 0.0, 0.0 });
        return mesh;
    }

    [[nodiscard]] TriangleMesh MakeDegenerateSection()
    {
        TriangleMesh mesh;
        const IndexType first = mesh.AddVertex({ 1.0, 0.0, 0.0 });
        const IndexType second = mesh.AddVertex({ 1.0, 1.0, 0.0 });
        const IndexType third = mesh.AddVertex({ 1.0, 2.0, 0.0 });
        mesh.AddTriangle({ first, second, third }, { 1.0, 0.0, 0.0 });
        return mesh;
    }

}  // namespace

TEST_SUITE("LocalSurfaceClassifierTestSuite")
{

    TEST_CASE("LocalSurfaceClassifierTestSuite - Classifies Both Sides Independently Of Section Scale")
    {
        for (const double scale : std::array{ 1.0, 1e-6 }) {
            const TriangleMesh section = MakeOrientedSection(scale);
            const GeometryTolerance tolerance =
                GeometryTolerance::FromScale({ .length_scale = scale, .coordinate_scale = 1.0 });
            const PointType inside{ 0.0, scale / 3.0, scale / 3.0 };
            const PointType outside{ 2.0 * scale, scale / 3.0, scale / 3.0 };
            CHECK(
                embedding::detail::ClassifyOnBoundedSide(inside, section.View(), tolerance.ZeroLength())
                == embedding::detail::LocalSurfaceClassification::Inside
            );
            CHECK(
                embedding::detail::ClassifyOnBoundedSide(outside, section.View(), tolerance.ZeroLength())
                == embedding::detail::LocalSurfaceClassification::Outside
            );
        }
    }

    TEST_CASE("LocalSurfaceClassifierTestSuite - Inconclusive Evidence Remains Inconclusive")
    {
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const TriangleMesh empty_section;
        CHECK(
            embedding::detail::ClassifyOnBoundedSide(PointType{}, empty_section.View(), tolerance.ZeroLength())
            == embedding::detail::LocalSurfaceClassification::Inconclusive
        );

        const TriangleMesh degenerate_section = MakeDegenerateSection();
        CHECK(
            embedding::detail::ClassifyOnBoundedSide(
                PointType{ 0.0, 0.5, 0.0 }, degenerate_section.View(), tolerance.ZeroLength()
            )
            == embedding::detail::LocalSurfaceClassification::Inconclusive
        );
    }

    TEST_CASE("LocalSurfaceClassifierTestSuite - Tri State Classification Is Repeatable And Concurrent")
    {
        const TriangleMesh section = MakeOrientedSection(1.0);
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const PointType point{ 0.0, 1.0 / 3.0, 1.0 / 3.0 };
        int failures = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : failures)
#endif
        for (int i = 0; i < 100; ++i) {
            if (embedding::detail::ClassifyOnBoundedSide(point, section.View(), tolerance.ZeroLength())
                != embedding::detail::LocalSurfaceClassification::Inside) {
                ++failures;
            }
        }
        CHECK((failures) == (0));
    }

}  // TEST_SUITE

}  // namespace queso::testing

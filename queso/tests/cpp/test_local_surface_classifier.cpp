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
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/embedding/local_surface_classifier.h"
#include "queso/includes/checks.hpp"

namespace queso::Testing {
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

BOOST_AUTO_TEST_SUITE(LocalSurfaceClassifierTestSuite)

BOOST_AUTO_TEST_CASE(ClassifiesBothSidesIndependentlyOfSectionScale)
{
    for (const double scale : std::array{ 1.0, 1e-6 }) {
        const TriangleMesh section = MakeOrientedSection(scale);
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = scale, .coordinate_scale = 1.0 });
        const PointType inside{ 0.0, scale / 3.0, scale / 3.0 };
        const PointType outside{ 2.0 * scale, scale / 3.0, scale / 3.0 };
        QuESo_CHECK(
            embedding::detail::ClassifyOnBoundedSide(inside, section.View(), tolerance.ZeroLength())
            == embedding::detail::LocalSurfaceClassification::Inside
        );
        QuESo_CHECK(
            embedding::detail::ClassifyOnBoundedSide(outside, section.View(), tolerance.ZeroLength())
            == embedding::detail::LocalSurfaceClassification::Outside
        );
    }
}

BOOST_AUTO_TEST_CASE(InconclusiveEvidenceRemainsInconclusive)
{
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const TriangleMesh empty_section;
    QuESo_CHECK(
        embedding::detail::ClassifyOnBoundedSide(PointType{}, empty_section.View(), tolerance.ZeroLength())
        == embedding::detail::LocalSurfaceClassification::Inconclusive
    );

    const TriangleMesh degenerate_section = MakeDegenerateSection();
    QuESo_CHECK(
        embedding::detail::ClassifyOnBoundedSide(
            PointType{ 0.0, 0.5, 0.0 }, degenerate_section.View(), tolerance.ZeroLength()
        )
        == embedding::detail::LocalSurfaceClassification::Inconclusive
    );
}

BOOST_AUTO_TEST_CASE(TriStateClassificationIsRepeatableAndConcurrent)
{
    const TriangleMesh section = MakeOrientedSection(1.0);
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
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
    QuESo_CHECK_EQUAL(failures, 0);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

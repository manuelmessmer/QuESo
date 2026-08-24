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
#include "queso/embedding/cell_surface_section.h"
#include "queso/includes/checks.hpp"

namespace queso::Testing {

namespace {

    [[nodiscard]] TriangleMesh MakeSurface()
    {
        TriangleMesh result;
        result.AddVertex({ 0.0, 0.0, 0.0 });
        result.AddVertex({ 1.0, 0.0, 0.0 });
        result.AddVertex({ 0.0, 1.0, 0.0 });
        result.AddTriangle({ 0, 1, 2 }, { 0.0, 0.0, 1.0 });
        return result;
    }

    [[nodiscard]] embedding::CellFaceContours MakeContours()
    {
        embedding::CellFaceContours result;
        result.segments[embedding::CellFace{ 0, false }.Index()].push_back(
            { { 0.0, 0.25, 0.0 }, { 0.0, 0.75, 0.0 }, { 0.0, 0.0, 1.0 }, 7 }
        );
        return result;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(CellSurfaceSectionTestSuite)

BOOST_AUTO_TEST_CASE(CellFaceIndexingIsDenseAndReversible)
{
    for (IndexType index = 0; index < 6; ++index) {
        const embedding::CellFace face = embedding::CellFace::FromIndex(index);
        QuESo_CHECK_EQUAL(face.Index(), index);
        QuESo_CHECK_EQUAL(face.axis, index / 2);
        QuESo_CHECK(face.is_upper == (index % 2 == 1));
    }
}

BOOST_AUTO_TEST_CASE(OwnsSurfaceContoursBoundsAndTolerance)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    embedding::CellSurfaceSection section(MakeSurface(), MakeContours(), bounds, tolerance);

    QuESo_CHECK_EQUAL(section.SurfaceView().NumOfTriangles(), 1);
    QuESo_CHECK_EQUAL(section.SurfaceView().NumOfTriangles(), 1);
    QuESo_CHECK_EQUAL(section.FaceSegments({ 0, false }).size(), 1);
    QuESo_CHECK_EQUAL(section.FaceSegments({ 0, false })[0].source_triangle, 7);
    QuESo_CHECK(section.FaceSegments({ 0, true }).empty());
    QuESo_CHECK_POINT_NEAR(section.CellBounds().lower, bounds.lower, 0.0);
    QuESo_CHECK_POINT_NEAR(section.CellBounds().upper, bounds.upper, 0.0);
    QuESo_CHECK_EQUAL(section.GetGeometryTolerance().SnapDistance(), tolerance.SnapDistance());

    embedding::CellSurfaceSection moved(std::move(section));
    QuESo_CHECK_EQUAL(moved.SurfaceView().NumOfTriangles(), 1);
    QuESo_CHECK_EQUAL(moved.FaceSegments({ 0, false }).size(), 1);
}

BOOST_AUTO_TEST_CASE(RejectsInvalidBoundsAndToleranceSnapshot)
{
    if constexpr (NOTDEBUG) { return; }
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    BOOST_CHECK_THROW(
        embedding::CellSurfaceSection(
            MakeSurface(), MakeContours(), MakeBox({ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 1.0 }), tolerance
        ),
        Exception
    );

    const GeometryTolerance large_tolerance =
        GeometryTolerance::FromScale({ .length_scale = 1e6, .coordinate_scale = 1e6 });
    BOOST_CHECK_THROW(
        embedding::CellSurfaceSection(
            MakeSurface(), MakeContours(), MakeBox({ 0.0, 0.0, 0.0 }, { 1e-7, 1.0, 1.0 }), large_tolerance
        ),
        Exception
    );
}

BOOST_AUTO_TEST_CASE(RejectsInvalidSegmentGeometry)
{
    if constexpr (NOTDEBUG) { return; }
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

    auto off_plane = MakeContours();
    off_plane.segments[0][0].first[0] = tolerance.SnapDistance();
    BOOST_CHECK_THROW(embedding::CellSurfaceSection(MakeSurface(), std::move(off_plane), bounds, tolerance), Exception);

    auto outside = MakeContours();
    outside.segments[0][0].second[1] = 2.0;
    BOOST_CHECK_THROW(embedding::CellSurfaceSection(MakeSurface(), std::move(outside), bounds, tolerance), Exception);

    auto non_finite = MakeContours();
    non_finite.segments[0][0].source_normal[0] = std::numeric_limits<double>::quiet_NaN();
    BOOST_CHECK_THROW(
        embedding::CellSurfaceSection(MakeSurface(), std::move(non_finite), bounds, tolerance), Exception
    );
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

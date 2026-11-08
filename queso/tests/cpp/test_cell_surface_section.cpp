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
#include "queso/embedding/cell_surface_section.h"

namespace queso::testing {

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

TEST_SUITE("CellSurfaceSectionTestSuite")
{

    TEST_CASE("CellSurfaceSectionTestSuite - Cell Face Indexing Is Dense And Reversible")
    {
        for (IndexType index = 0; index < 6; ++index) {
            const embedding::CellFace face = embedding::CellFace::FromIndex(index);
            CHECK((face.Index()) == (index));
            CHECK((face.axis) == (index / 2));
            CHECK(face.is_upper == (index % 2 == 1));
        }
    }

    TEST_CASE("CellSurfaceSectionTestSuite - Owns Surface Contours Bounds And Tolerance")
    {
        const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        embedding::CellSurfaceSection section(MakeSurface(), MakeContours(), bounds, tolerance);

        CHECK((section.SurfaceView().NumOfTriangles()) == (1));
        CHECK((section.SurfaceView().NumOfTriangles()) == (1));
        REQUIRE((section.FaceSegments({ 0, false }).size()) == (1));
        CHECK((section.FaceSegments({ 0, false })[0].source_triangle) == (7));
        REQUIRE(section.FaceSegments({ 0, true }).empty());
        CheckPointNear(section.CellBounds().lower, bounds.lower, 0.0);
        CheckPointNear(section.CellBounds().upper, bounds.upper, 0.0);
        CHECK((section.GetGeometryTolerance().SnapDistance()) == (tolerance.SnapDistance()));

        embedding::CellSurfaceSection moved(std::move(section));
        CHECK((moved.SurfaceView().NumOfTriangles()) == (1));
        REQUIRE((moved.FaceSegments({ 0, false }).size()) == (1));
    }

    TEST_CASE("CellSurfaceSectionTestSuite - Rejects Invalid Bounds And Tolerance Snapshot")
    {
        if (!QuESo_DIAGNOSTICS_ENABLED()) { return; }
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        CHECK_THROWS_AS(
            embedding::CellSurfaceSection(
                MakeSurface(), MakeContours(), MakeBox({ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 1.0 }), tolerance
            ),
            Exception
        );

        const GeometryTolerance large_tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1e6, .coordinate_scale = 1e6 });
        CHECK_THROWS_AS(
            embedding::CellSurfaceSection(
                MakeSurface(), MakeContours(), MakeBox({ 0.0, 0.0, 0.0 }, { 1e-7, 1.0, 1.0 }), large_tolerance
            ),
            Exception
        );
    }

    TEST_CASE("CellSurfaceSectionTestSuite - Rejects Invalid Segment Geometry")
    {
        if (!QuESo_DIAGNOSTICS_ENABLED()) { return; }
        const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

        auto off_plane = MakeContours();
        off_plane.segments[0][0].first[0] = tolerance.SnapDistance();
        CHECK_THROWS_AS(
            embedding::CellSurfaceSection(MakeSurface(), std::move(off_plane), bounds, tolerance), Exception
        );

        auto outside = MakeContours();
        outside.segments[0][0].second[1] = 2.0;
        CHECK_THROWS_AS(embedding::CellSurfaceSection(MakeSurface(), std::move(outside), bounds, tolerance), Exception);

        auto non_finite = MakeContours();
        non_finite.segments[0][0].source_normal[0] = std::numeric_limits<double>::quiet_NaN();
        CHECK_THROWS_AS(
            embedding::CellSurfaceSection(MakeSurface(), std::move(non_finite), bounds, tolerance), Exception
        );
    }

    TEST_CASE("CellSurfaceSectionTestSuite - Rejects Non Finite Surface Data When Diagnostics Are Enabled")
    {
        if (!QuESo_DIAGNOSTICS_ENABLED()) { return; }
        const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
        const double nan = std::numeric_limits<double>::quiet_NaN();
        TriangleMesh invalid_vertices;
        invalid_vertices.AddVertex({ nan, 0.0, 0.0 });
        CHECK_THROWS_AS(embedding::CellSurfaceSection(std::move(invalid_vertices), {}, bounds, tolerance), Exception);
        TriangleMesh invalid_normals;
        invalid_normals.AddVertex({ 0.0, 0.0, 0.0 });
        invalid_normals.AddVertex({ 1.0, 0.0, 0.0 });
        invalid_normals.AddVertex({ 0.0, 1.0, 0.0 });
        invalid_normals.AddTriangle({ 0, 1, 2 }, { nan, 0.0, 1.0 });
        CHECK_THROWS_AS(embedding::CellSurfaceSection(std::move(invalid_normals), {}, bounds, tolerance), Exception);
    }

}  // TEST_SUITE

}  // namespace queso::testing

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
#include <algorithm>
#include <array>
#include <vector>

//// Project includes
#include "queso/embedding/cell_face_closure.h"
#include "queso/includes/checks.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {
namespace {

    [[nodiscard]] PointType Lift(embedding::CellFace Face, const std::array<double, 2>& rPoint)
    {
        PointType result{};
        result[Face.axis] = Face.is_upper ? 1.0 : 0.0;
        result[(Face.axis + 1) % 3] = rPoint[0];
        result[(Face.axis + 2) % 3] = rPoint[1];
        return result;
    }

    [[nodiscard]] embedding::CellFaceSegment MakeSegment(
        embedding::CellFace Face,
        const std::array<double, 2>& rFirst,
        const std::array<double, 2>& rSecond,
        IndexType SourceTriangle
    )
    {
        const PointType first = Lift(Face, rFirst);
        const PointType second = Lift(Face, rSecond);
        PointType outward{};
        outward[Face.axis] = Face.is_upper ? 1.0 : -1.0;
        Vector3d source_normal = Math::Cross(second - first, outward);
        source_normal /= Math::Norm(source_normal);
        return { first, second, source_normal, SourceTriangle };
    }

    [[nodiscard]] std::vector<embedding::CellFaceSegment>
        MakeLoop(embedding::CellFace Face, std::span<const std::array<double, 2>> rPoints)
    {
        std::vector<embedding::CellFaceSegment> result;
        for (IndexType i = 0; i < rPoints.size(); ++i) {
            const IndexType first = Face.is_upper ? i : rPoints.size() - 1 - i;
            const IndexType second =
                Face.is_upper ? (i + 1) % rPoints.size() : (rPoints.size() - 2 - i + rPoints.size()) % rPoints.size();
            result.push_back(MakeSegment(Face, rPoints[first], rPoints[second], i));
        }
        return result;
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(CellFaceClosureTestSuite)

BOOST_AUTO_TEST_CASE(FillsOrEmptiesUncutFaces)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    for (IndexType face_index = 0; face_index < 6; ++face_index) {
        const embedding::CellFace face = embedding::CellFace::FromIndex(face_index);
        const TriangleMesh empty =
            embedding::detail::BuildCellFaceClosure({}, bounds, face, tolerance, [](PointView) { return false; });
        QuESo_CHECK_EQUAL(empty.NumOfTriangles(), IndexType{ 0 });

        const TriangleMesh full =
            embedding::detail::BuildCellFaceClosure({}, bounds, face, tolerance, [](PointView) { return true; });
        QuESo_CHECK_EQUAL(full.NumOfTriangles(), IndexType{ 2 });
        QuESo_CHECK_NEAR(MeshUtilities::Area(full.View()), 1.0, EPS4);
        PointType outward{};
        outward[face.axis] = face.is_upper ? 1.0 : -1.0;
        for (const auto& r_triangle : full.Triangles<WithNormals>()) {
            QuESo_CHECK_POINT_NEAR(r_triangle.Normal, outward, 0.0);
            QuESo_CHECK_GT(
                Math::Dot(Math::Cross(r_triangle.P2 - r_triangle.P1, r_triangle.P3 - r_triangle.P1), outward), 0.0
            );
        }
    }
}

BOOST_AUTO_TEST_CASE(TriangulatesConvexLoopOnAllFaces)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    constexpr std::array points{
        std::array{ 0.2, 0.2 }, std::array{ 0.8, 0.2 }, std::array{ 0.8, 0.8 }, std::array{ 0.2, 0.8 }
    };
    for (IndexType face_index = 0; face_index < 6; ++face_index) {
        const embedding::CellFace face = embedding::CellFace::FromIndex(face_index);
        const auto segments = MakeLoop(face, points);
        const TriangleMesh closure =
            embedding::detail::BuildCellFaceClosure(segments, bounds, face, tolerance, [](PointView) { return false; });
        QuESo_CHECK_NEAR(MeshUtilities::Area(closure.View()), 0.36, EPS4);
    }
}

BOOST_AUTO_TEST_CASE(TriangulatesOppositeAndAdjacentBoundaryChords)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const embedding::CellFace face{ 2, true };

    const std::array horizontal{ MakeSegment(face, { 0.0, 0.5 }, { 1.0, 0.5 }, 0) };
    const TriangleMesh upper =
        embedding::detail::BuildCellFaceClosure(horizontal, bounds, face, tolerance, [](PointView rPoint) {
            return rPoint[1] > 0.5;
        });
    QuESo_CHECK_NEAR(MeshUtilities::Area(upper.View()), 0.5, EPS4);

    const std::array adjacent{ MakeSegment(face, { 0.0, 0.5 }, { 0.5, 0.0 }, 0) };
    const TriangleMesh large_side =
        embedding::detail::BuildCellFaceClosure(adjacent, bounds, face, tolerance, [](PointView rPoint) {
            return rPoint[0] + rPoint[1] > 0.5;
        });
    QuESo_CHECK_NEAR(MeshUtilities::Area(large_side.View()), 0.875, EPS4);
}

BOOST_AUTO_TEST_CASE(TriangulatesCavityAndDisconnectedLoops)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const embedding::CellFace face{ 2, true };
    constexpr std::array outer{
        std::array{ 0.1, 0.1 }, std::array{ 0.9, 0.1 }, std::array{ 0.9, 0.9 }, std::array{ 0.1, 0.9 }
    };
    constexpr std::array hole{
        std::array{ 0.3, 0.3 }, std::array{ 0.3, 0.7 }, std::array{ 0.7, 0.7 }, std::array{ 0.7, 0.3 }
    };
    auto segments = MakeLoop(face, outer);
    auto hole_segments = MakeLoop(face, hole);
    segments.insert(segments.end(), hole_segments.begin(), hole_segments.end());
    const TriangleMesh cavity =
        embedding::detail::BuildCellFaceClosure(segments, bounds, face, tolerance, [](PointView rPoint) {
            const bool in_outer = rPoint[0] > 0.1 && rPoint[0] < 0.9 && rPoint[1] > 0.1 && rPoint[1] < 0.9;
            const bool in_hole = rPoint[0] > 0.3 && rPoint[0] < 0.7 && rPoint[1] > 0.3 && rPoint[1] < 0.7;
            return in_outer && !in_hole;
        });
    QuESo_CHECK_NEAR(MeshUtilities::Area(cavity.View()), 0.48, EPS4);
}

BOOST_AUTO_TEST_CASE(TriangulatesConcaveLoopAndIsIdempotent)
{
    const BoundingBoxType bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const embedding::CellFace face{ 2, true };
    constexpr std::array points{ std::array{ 0.1, 0.1 }, std::array{ 0.9, 0.1 }, std::array{ 0.9, 0.4 },
                                 std::array{ 0.4, 0.4 }, std::array{ 0.4, 0.9 }, std::array{ 0.1, 0.9 } };
    const auto segments = MakeLoop(face, points);
    const auto classifier = [](PointView rPoint) {
        return (rPoint[0] > 0.1 && rPoint[0] < 0.9 && rPoint[1] > 0.1 && rPoint[1] < 0.4)
               || (rPoint[0] > 0.1 && rPoint[0] < 0.4 && rPoint[1] >= 0.4 && rPoint[1] < 0.9);
    };
    const TriangleMesh first = embedding::detail::BuildCellFaceClosure(segments, bounds, face, tolerance, classifier);
    const TriangleMesh second = embedding::detail::BuildCellFaceClosure(segments, bounds, face, tolerance, classifier);
    QuESo_CHECK_NEAR(MeshUtilities::Area(first.View()), 0.39, EPS4);
    QuESo_CHECK(std::ranges::equal(first.Vertices(), second.Vertices()));
    QuESo_CHECK(std::ranges::equal(first.TriangleIndices(), second.TriangleIndices()));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

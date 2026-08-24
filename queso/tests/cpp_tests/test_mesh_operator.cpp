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
#include <future>
#include <memory>

//// External includes
#include <boost/test/unit_test.hpp>

//// Project includes
#include "queso/embedding/mesh_operator.h"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::Testing {
BOOST_AUTO_TEST_SUITE(MeshOperatorTestSuite)

BOOST_AUTO_TEST_CASE(OperatorOwnsOneStableQuery)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshOperator mesh_operator(
        mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
    );
    const embedding::MeshQuery* p_query = &mesh_operator.Query();

    QuESo_CHECK_EQUAL(p_query, &mesh_operator.Query());
    QuESo_CHECK(p_query->IsWithinBoundingBox(PointType{ 0.5, 0.5, 0.5 }));
    QuESo_CHECK_EQUAL(p_query->MeshView().NumOfTriangles(), mesh.NumOfTriangles());
}

BOOST_AUTO_TEST_CASE(PointClassificationPreservesBoxBehavior)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshOperator mesh_operator(
        mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
    );

    QuESo_CHECK(mesh_operator.IsInside(PointType{ 0.5, 0.5, 0.5 }));
    QuESo_CHECK(!mesh_operator.IsInside(PointType{ 2.0, 2.0, 2.0 }));
}

BOOST_AUTO_TEST_CASE(PointClassificationIsDeterministicAndThreadSafe)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const auto p_mesh_operator = std::make_shared<const embedding::MeshOperator>(
        mesh.View(), GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 })
    );
    constexpr PointType inside{ 0.37, 0.41, 0.53 };
    constexpr PointType outside{ 2.0, 2.0, 2.0 };
    constexpr IndexType number_of_workers = 8;
    std::array<std::future<bool>, number_of_workers> results;
    for (auto& r_result : results) {
        r_result = std::async(std::launch::async, [p_mesh_operator, inside, outside]() {
            for (IndexType repetition = 0; repetition < 100; ++repetition) {
                if (!p_mesh_operator->IsInside(inside) || p_mesh_operator->IsInside(outside)) { return false; }
            }
            return true;
        });
    }
    for (auto& r_result : results) { QuESo_CHECK(r_result.get()); }
}

BOOST_AUTO_TEST_CASE(ExactBoundDomainClippingPreservesProvenance)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const BoundingBoxType bounds = MakeBox({ -0.1, 0.2, 0.2 }, { 0.6, 0.8, 0.8 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 0.7, .coordinate_scale = 1.0 });
    const embedding::MeshOperator mesh_operator(mesh.View(), tolerance);
    const auto candidates = mesh_operator.Query().GetAabbCandidates(bounds.lower, bounds.upper);

    const auto section = mesh_operator.ClipCellSurfaceSection(candidates, bounds);
    QuESo_CHECK(section.SurfaceView().NumOfTriangles() > 0);
    bool has_segment = false;
    for (IndexType face = 0; face < 6; ++face) {
        for (const auto& r_segment : section.FaceSegments(embedding::CellFace::FromIndex(face))) {
            has_segment = true;
            QuESo_CHECK(r_segment.source_triangle < mesh.NumOfTriangles());
            for (const PointType& r_point : { r_segment.first, r_segment.second }) {
                for (IndexType axis = 0; axis < 3; ++axis) {
                    QuESo_CHECK(r_point[axis] >= bounds.lower[axis]);
                    QuESo_CHECK(r_point[axis] <= bounds.upper[axis]);
                }
            }
        }
    }
    QuESo_CHECK(has_segment);
}

BOOST_AUTO_TEST_CASE(WiderCandidateDiscoveryKeepsExactClippingBounds)
{
    TriangleMesh mesh;
    const IndexType first = mesh.AddVertex({ -1.0, 0.2, 0.5 });
    const IndexType second = mesh.AddVertex({ 1.0, 0.2, 0.5 });
    const IndexType third = mesh.AddVertex({ 1.0, 0.8, 0.5 });
    mesh.AddTriangle({ first, second, third }, { 0.0, 0.0, 1.0 });
    const BoundingBoxType exact_bounds = MakeBox({ 0.0, 0.0, 0.0 }, { 0.5, 1.0, 1.0 });
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });
    const embedding::MeshOperator mesh_operator(mesh.View(), tolerance);
    constexpr std::array<IndexType, 0> no_candidates{};
    constexpr std::array<IndexType, 1> wider_candidates{ 0 };

    const auto empty = mesh_operator.ClipCellSurfaceSection(no_candidates, exact_bounds);
    QuESo_CHECK_EQUAL(empty.SurfaceView().NumOfTriangles(), IndexType{ 0 });
    const auto recovered = mesh_operator.ClipCellSurfaceSection(wider_candidates, exact_bounds);
    QuESo_CHECK(recovered.SurfaceView().NumOfTriangles() > 0);
    const auto CheckVertex = [&](PointView rVertex) {
        for (IndexType axis = 0; axis < 3; ++axis) {
            QuESo_CHECK(rVertex[axis] >= exact_bounds.lower[axis]);
            QuESo_CHECK(rVertex[axis] <= exact_bounds.upper[axis]);
        }
    };
    recovered.SurfaceView().VisitEachTriangle<WithoutNormals>([&](const auto& rTriangle) {
        CheckVertex(rTriangle.P1);
        CheckVertex(rTriangle.P2);
        CheckVertex(rTriangle.P3);
    });
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

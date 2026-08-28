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
#include <limits>
#include <vector>

//// Project includes
#include "queso/embedding/aabb_triangle_intersector.h"
#include "queso/embedding/mesh_query.h"
#include "queso/includes/checks.hpp"
#include "queso/utilities/mesh_utilities.h"

namespace queso::Testing {
namespace {

    using embedding::Ray;

    [[nodiscard]] GeometryTolerance MakeTolerance(double LengthScale = 1.0)
    { return GeometryTolerance::FromScale({ .length_scale = LengthScale, .coordinate_scale = 1.0 }); }

    void AddTriangle(TriangleMesh& rMesh, PointView rFirst, PointView rSecond, PointView rThird)
    {
        const IndexType first = rMesh.AddVertex({ rFirst[0], rFirst[1], rFirst[2] });
        const IndexType second = rMesh.AddVertex({ rSecond[0], rSecond[1], rSecond[2] });
        const IndexType third = rMesh.AddVertex({ rThird[0], rThird[1], rThird[2] });
        rMesh.AddTriangle({ first, second, third }, { 0.0, 0.0, 1.0 });
    }

}  // namespace

BOOST_AUTO_TEST_SUITE(MeshQueryTestSuite)

BOOST_AUTO_TEST_CASE(CandidatesAreConservativeAndExactTestsReuseThem)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const PointType lower{ -0.1, 0.2, 0.2 };
    const PointType upper{ 0.1, 0.8, 0.8 };

    const auto candidates = query.GetAabbCandidates(lower, upper);
    const auto intersected = query.GetIntersectedTriangleIds(candidates, lower, upper);
    const embedding::detail::AabbTriangleIntersector aabb(lower, upper);

    mesh.View().VisitEachTriangle<WithoutNormals>([&, triangle_id = IndexType{ 0 }](const auto& rTriangle) mutable {
        if (aabb.IntersectsTriangle(rTriangle, 0.0)) {
            QuESo_CHECK(std::ranges::find(candidates, triangle_id) != candidates.end());
        }
        ++triangle_id;
    });

    for (const IndexType triangle_id : intersected) {
        QuESo_CHECK(std::ranges::find(candidates, triangle_id) != candidates.end());
    }
    QuESo_CHECK(!intersected.empty());
    QuESo_CHECK(query.IntersectsAabb(candidates, lower, upper));
}

BOOST_AUTO_TEST_CASE(EmptyCandidatesDoNotIntersect)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const PointType lower{ -0.1, 0.2, 0.2 };
    const PointType upper{ 0.1, 0.8, 0.8 };
    const std::array<IndexType, 0> candidates{};

    QuESo_CHECK(!query.IntersectsAabb(candidates, lower, upper));
    QuESo_CHECK(query.GetIntersectedTriangleIds(candidates, lower, upper).empty());
}

BOOST_AUTO_TEST_CASE(TouchingRespectsExactTolerance)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const PointType lower{ 1.0, 0.2, 0.2 };
    const PointType upper{ 2.0, 0.8, 0.8 };
    const auto candidates = query.GetAabbCandidates(lower, upper);

    QuESo_CHECK(query.IntersectsAabb(candidates, lower, upper));
    QuESo_CHECK(!query.IntersectsAabb(candidates, lower, upper, AabbIntersectionPolicy::SnapEroded));
    QuESo_CHECK(!query.GetIntersectedTriangleIds(candidates, lower, upper).empty());
}

BOOST_AUTO_TEST_CASE(ReportsMeshBoundsAndClosedRayClassification)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ -1.0, -2.0, -3.0 }, { 2.0, 3.0, 4.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(7.0));
    const auto bounds = query.BoundingBox();
    BOOST_REQUIRE(bounds.has_value());

    QuESo_CHECK_EQUAL(bounds->lower, (PointType{ -1.0, -2.0, -3.0 }));
    QuESo_CHECK_EQUAL(bounds->upper, (PointType{ 2.0, 3.0, 4.0 }));

    const Ray ray(PointType{ 0.0, 0.0, 0.0 }, Vector3d{ 1.0, 0.123, 0.312 });
    const auto [is_inside, success] = query.Classify(ray);
    QuESo_CHECK(success);
    QuESo_CHECK(is_inside);
}

BOOST_AUTO_TEST_CASE(OpenRayClassificationUsesTriangleOrientation)
{
    TriangleMesh mesh;
    mesh.AddVertex({ 0.0, 0.0, 0.0 });
    mesh.AddVertex({ 1.0, 0.0, 0.0 });
    mesh.AddVertex({ 0.0, 1.0, 0.0 });
    mesh.AddTriangle({ 0, 1, 2 }, { 0.0, 0.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::OrientedSurface, MakeTolerance());

    const auto [is_inside, inside_success] =
        query.Classify(Ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
    const auto [is_outside, outside_success] =
        query.Classify(Ray(PointType{ 0.2, 0.2, 1.0 }, Vector3d{ 0.0, 0.0, -1.0 }));

    QuESo_CHECK(inside_success);
    QuESo_CHECK(outside_success);
    QuESo_CHECK(is_inside);
    QuESo_CHECK(!is_outside);
}

BOOST_AUTO_TEST_CASE(ExactFilteringRejectsFalsePositivesAndPreservesCandidateOrder)
{
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 2.0, 0.0, 0.0 }, PointType{ 0.0, 2.0, 0.0 });
    AddTriangle(mesh, PointType{ 1.55, 1.55, 0.0 }, PointType{ 1.75, 1.55, 0.0 }, PointType{ 1.55, 1.75, 0.0 });
    AddTriangle(mesh, PointType{ 1.6, 1.6, 0.05 }, PointType{ 1.7, 1.6, 0.05 }, PointType{ 1.6, 1.7, 0.05 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(2.0));
    constexpr PointType lower{ 1.5, 1.5, -0.1 };
    constexpr PointType upper{ 1.8, 1.8, 0.1 };
    const auto conservative = query.GetAabbCandidates(lower, upper);
    QuESo_CHECK(std::ranges::find(conservative, IndexType{ 0 }) != conservative.end());
    constexpr std::array<IndexType, 1> false_positive{ 0 };
    QuESo_CHECK(!query.IntersectsAabb(false_positive, lower, upper));
    constexpr std::array<IndexType, 3> candidates{ 2, 0, 1 };
    const auto exact = query.GetIntersectedTriangleIds(candidates, lower, upper);
    QuESo_CHECK(exact == (std::vector<IndexType>{ 2, 1 }));
}

BOOST_AUTO_TEST_CASE(ClosedClassificationCoversOutsideAndAmbiguousRays)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const auto [outside, outside_success] =
        query.Classify(Ray(PointType{ -1.0, 0.23, 0.31 }, Vector3d{ 1.0, 0.01, 0.02 }));
    QuESo_CHECK(outside_success);
    QuESo_CHECK(!outside);
    const auto [miss, miss_success] = query.Classify(Ray(PointType{ 2.0, 2.0, 2.0 }, Vector3d{ 1.0, 0.1, 0.2 }));
    QuESo_CHECK(miss_success);
    QuESo_CHECK(!miss);
    const auto [edge, edge_success] = query.Classify(Ray(PointType{ 0.5, 0.5, 0.5 }, Vector3d{ 1.0, 0.2, 0.2 }));
    QuESo_CHECK(!edge);
    QuESo_CHECK(!edge_success);
}

BOOST_AUTO_TEST_CASE(ParallelClosedCandidateIsInconclusiveWithoutReadingHitCoordinates)
{
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const auto [is_inside, is_conclusive] = query.Classify(Ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.1, 0.0 }));
    QuESo_CHECK(!is_inside);
    QuESo_CHECK(!is_conclusive);
}

BOOST_AUTO_TEST_CASE(BoundaryOriginsClassifyOutsideConclusive)
{
    const TriangleMesh box = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery closed_query(box.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    const auto [closed_inside, closed_success] =
        closed_query.Classify(Ray(PointType{ 0.0, 0.23, 0.31 }, Vector3d{ 1.0, 0.01, 0.02 }));
    QuESo_CHECK(!closed_inside);
    QuESo_CHECK(closed_success);

    TriangleMesh surface;
    AddTriangle(surface, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
    const embedding::MeshQuery oriented_query(
        surface.View(), embedding::MeshQueryMode::OrientedSurface, MakeTolerance()
    );
    const auto [surface_inside, surface_success] =
        oriented_query.Classify(Ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
    QuESo_CHECK(!surface_inside);
    QuESo_CHECK(surface_success);
}

BOOST_AUTO_TEST_CASE(OrientedSurfaceUsesNearestHitAndDistinguishesMissFromEdge)
{
    TriangleMesh mesh;
    AddTriangle(mesh, PointType{ 0.0, 1.0, 1.0 }, PointType{ 1.0, 0.0, 1.0 }, PointType{ 0.0, 0.0, 1.0 });
    AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::OrientedSurface, MakeTolerance());
    const auto [inside, inside_success] = query.Classify(Ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
    QuESo_CHECK(inside_success);
    QuESo_CHECK(inside);
    const auto [edge, edge_success] = query.Classify(Ray(PointType{ 0.5, 0.5, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
    QuESo_CHECK(!edge);
    QuESo_CHECK(!edge_success);
    const auto [miss, miss_success] = query.Classify(Ray(PointType{ 2.0, 2.0, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
    QuESo_CHECK(!miss);
    QuESo_CHECK(miss_success);
}

BOOST_AUTO_TEST_CASE(ClosedClassificationSupportsSmallScaleMeshes)
{
    constexpr double scale = 1e-6;
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { scale, scale, scale });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(scale));
    const auto [is_inside, is_conclusive] =
        query.Classify(Ray(PointType{ 0.5 * scale, 0.5 * scale, 0.5 * scale }, Vector3d{ 1.0, 0.123, 0.312 }));
    QuESo_CHECK(is_conclusive);
    QuESo_CHECK(is_inside);
}

BOOST_AUTO_TEST_CASE(BoundingBoxAndInvalidInputContracts)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    QuESo_CHECK(query.IsWithinBoundingBox(PointType{ 0.0, 0.0, 0.0 }));
    QuESo_CHECK(query.IsWithinBoundingBox(PointType{ 1.0, 1.0, 1.0 }));
    QuESo_CHECK(!query.IsWithinBoundingBox(PointType{ -1e-12, 0.5, 0.5 }));
    QuESo_CHECK(!query.IsWithinBoundingBox(PointType{ std::numeric_limits<double>::quiet_NaN(), 0.5, 0.5 }));
    if constexpr (!NOTDEBUG) {
        BOOST_CHECK_THROW(
            (void)query.GetAabbCandidates(PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 1.0 }), Exception
        );
    }
    const auto candidates = query.GetAabbCandidates(PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 1.0, 1.0 });
    const TriangleMesh empty_mesh;
    const embedding::MeshQuery empty_query(empty_mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    QuESo_CHECK(!empty_query.BoundingBox().has_value());
    QuESo_CHECK(!empty_query.IsWithinBoundingBox(PointType{}));
}

BOOST_AUTO_TEST_CASE(RejectsInvalidRayDirections)
{
    BOOST_CHECK_THROW((void)Ray(PointType{}, Vector3d{}), queso::Exception);
    BOOST_CHECK_THROW(
        (void)Ray(PointType{}, Vector3d{ std::numeric_limits<double>::infinity(), 0.0, 0.0 }), queso::Exception
    );
}

BOOST_AUTO_TEST_CASE(DirectionMagnitudeDoesNotChangeClassification)
{
    const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
    const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
    constexpr PointType origin{ 0.5, 0.5, 0.5 };
    const auto unit_scale = query.Classify(Ray(origin, Vector3d{ 1.0, 0.123, 0.312 }));
    const auto large_scale = query.Classify(Ray(origin, Vector3d{ 10.0, 1.23, 3.12 }));
    QuESo_CHECK_EQUAL(unit_scale.is_inside, large_scale.is_inside);
    QuESo_CHECK_EQUAL(unit_scale.is_conclusive, large_scale.is_conclusive);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

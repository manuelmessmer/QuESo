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
#include <algorithm>
#include <array>
#include <limits>
#include <vector>

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/embedding/aabb_triangle_intersector.h"
#include "queso/embedding/mesh_query.h"
#include "queso/utilities/mesh_utilities.h"

namespace queso::testing {
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

TEST_SUITE("MeshQueryTestSuite")
{

    TEST_CASE("MeshQueryTestSuite - Candidates Are Conservative And Exact Tests Reuse Them")
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
                CHECK(std::ranges::find(candidates, triangle_id) != candidates.end());
            }
            ++triangle_id;
        });

        for (const IndexType triangle_id : intersected) {
            CHECK(std::ranges::find(candidates, triangle_id) != candidates.end());
        }
        REQUIRE(!intersected.empty());
        CHECK(query.IntersectsAabb(candidates, lower, upper));
    }

    TEST_CASE("MeshQueryTestSuite - Empty Candidates Do Not Intersect")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const PointType lower{ -0.1, 0.2, 0.2 };
        const PointType upper{ 0.1, 0.8, 0.8 };
        const std::array<IndexType, 0> candidates{};

        CHECK(!query.IntersectsAabb(candidates, lower, upper));
        REQUIRE(query.GetIntersectedTriangleIds(candidates, lower, upper).empty());
    }

    TEST_CASE("MeshQueryTestSuite - Touching Respects Exact Tolerance")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const PointType lower{ 1.0, 0.2, 0.2 };
        const PointType upper{ 2.0, 0.8, 0.8 };
        const auto candidates = query.GetAabbCandidates(lower, upper);

        CHECK(query.IntersectsAabb(candidates, lower, upper));
        CHECK(!query.IntersectsAabb(candidates, lower, upper, AabbIntersectionPolicy::SnapEroded));
        REQUIRE(!query.GetIntersectedTriangleIds(candidates, lower, upper).empty());
    }

    TEST_CASE("MeshQueryTestSuite - Reports Mesh Bounds And Closed Ray Classification")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ -1.0, -2.0, -3.0 }, { 2.0, 3.0, 4.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(7.0));
        const auto bounds = query.BoundingBox();
        REQUIRE(bounds.has_value());

        CHECK((bounds->lower) == ((PointType{ -1.0, -2.0, -3.0 })));
        CHECK((bounds->upper) == ((PointType{ 2.0, 3.0, 4.0 })));

        const Ray ray(PointType{ 0.0, 0.0, 0.0 }, Vector3d{ 1.0, 0.123, 0.312 });
        const auto [is_inside, success] = query.Classify(ray);
        CHECK(success);
        CHECK(is_inside);
    }

    TEST_CASE("MeshQueryTestSuite - Open Ray Classification Uses Triangle Orientation")
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

        CHECK(inside_success);
        CHECK(outside_success);
        CHECK(is_inside);
        CHECK(!is_outside);
    }

    TEST_CASE("MeshQueryTestSuite - Exact Filtering Rejects False Positives And Preserves Candidate Order")
    {
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 2.0, 0.0, 0.0 }, PointType{ 0.0, 2.0, 0.0 });
        AddTriangle(mesh, PointType{ 1.55, 1.55, 0.0 }, PointType{ 1.75, 1.55, 0.0 }, PointType{ 1.55, 1.75, 0.0 });
        AddTriangle(mesh, PointType{ 1.6, 1.6, 0.05 }, PointType{ 1.7, 1.6, 0.05 }, PointType{ 1.6, 1.7, 0.05 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(2.0));
        constexpr PointType lower{ 1.5, 1.5, -0.1 };
        constexpr PointType upper{ 1.8, 1.8, 0.1 };
        const auto conservative = query.GetAabbCandidates(lower, upper);
        CHECK(std::ranges::find(conservative, IndexType{ 0 }) != conservative.end());
        constexpr std::array<IndexType, 1> false_positive{ 0 };
        CHECK(!query.IntersectsAabb(false_positive, lower, upper));
        constexpr std::array<IndexType, 3> candidates{ 2, 0, 1 };
        const auto exact = query.GetIntersectedTriangleIds(candidates, lower, upper);
        CHECK(exact == (std::vector<IndexType>{ 2, 1 }));
    }

    TEST_CASE("MeshQueryTestSuite - Closed Classification Covers Outside And Ambiguous Rays")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const auto [outside, outside_success] =
            query.Classify(Ray(PointType{ -1.0, 0.23, 0.31 }, Vector3d{ 1.0, 0.01, 0.02 }));
        CHECK(outside_success);
        CHECK(!outside);
        const auto [miss, miss_success] = query.Classify(Ray(PointType{ 2.0, 2.0, 2.0 }, Vector3d{ 1.0, 0.1, 0.2 }));
        CHECK(miss_success);
        CHECK(!miss);
        const auto [edge, edge_success] = query.Classify(Ray(PointType{ 0.5, 0.5, 0.5 }, Vector3d{ 1.0, 0.2, 0.2 }));
        CHECK(!edge);
        CHECK(!edge_success);
    }

    TEST_CASE("MeshQueryTestSuite - Parallel Closed Candidate Is Inconclusive Without Reading Hit Coordinates")
    {
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const auto [is_inside, is_conclusive] =
            query.Classify(Ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.1, 0.0 }));
        CHECK(!is_inside);
        CHECK(!is_conclusive);
    }

    TEST_CASE("MeshQueryTestSuite - Boundary Origins Classify Outside Conclusive")
    {
        const TriangleMesh box = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery closed_query(box.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const auto [closed_inside, closed_success] =
            closed_query.Classify(Ray(PointType{ 0.0, 0.23, 0.31 }, Vector3d{ 1.0, 0.01, 0.02 }));
        CHECK(!closed_inside);
        CHECK(closed_success);

        TriangleMesh surface;
        AddTriangle(surface, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
        const embedding::MeshQuery oriented_query(
            surface.View(), embedding::MeshQueryMode::OrientedSurface, MakeTolerance()
        );
        const auto [surface_inside, surface_success] =
            oriented_query.Classify(Ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
        CHECK(!surface_inside);
        CHECK(surface_success);
    }

    TEST_CASE("MeshQueryTestSuite - Oriented Surface Uses Nearest Hit And Distinguishes Miss From Edge")
    {
        TriangleMesh mesh;
        AddTriangle(mesh, PointType{ 0.0, 1.0, 1.0 }, PointType{ 1.0, 0.0, 1.0 }, PointType{ 0.0, 0.0, 1.0 });
        AddTriangle(mesh, PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 0.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::OrientedSurface, MakeTolerance());
        const auto [inside, inside_success] =
            query.Classify(Ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
        CHECK(inside_success);
        CHECK(inside);
        const auto [edge, edge_success] = query.Classify(Ray(PointType{ 0.5, 0.5, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
        CHECK(!edge);
        CHECK(!edge_success);
        const auto [miss, miss_success] = query.Classify(Ray(PointType{ 2.0, 2.0, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 }));
        CHECK(!miss);
        CHECK(miss_success);
    }

    TEST_CASE("MeshQueryTestSuite - Closed Classification Supports Small Scale Meshes")
    {
        constexpr double scale = 1e-6;
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { scale, scale, scale });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance(scale));
        const auto [is_inside, is_conclusive] =
            query.Classify(Ray(PointType{ 0.5 * scale, 0.5 * scale, 0.5 * scale }, Vector3d{ 1.0, 0.123, 0.312 }));
        CHECK(is_conclusive);
        CHECK(is_inside);
    }

    TEST_CASE("MeshQueryTestSuite - Bounding Box And Invalid Input Diagnostics")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        CHECK(query.IsWithinBoundingBox(PointType{ 0.0, 0.0, 0.0 }));
        CHECK(query.IsWithinBoundingBox(PointType{ 1.0, 1.0, 1.0 }));
        CHECK(!query.IsWithinBoundingBox(PointType{ -1e-12, 0.5, 0.5 }));
        CHECK(!query.IsWithinBoundingBox(PointType{ std::numeric_limits<double>::quiet_NaN(), 0.5, 0.5 }));
        if (QuESo_DIAGNOSTICS_ENABLED()) {
            CHECK_THROWS_AS(
                (void)query.GetAabbCandidates(PointType{ 1.0, 0.0, 0.0 }, PointType{ 0.0, 1.0, 1.0 }), Exception
            );
        }
        const auto candidates = query.GetAabbCandidates(PointType{ 0.0, 0.0, 0.0 }, PointType{ 1.0, 1.0, 1.0 });
        const TriangleMesh empty_mesh;
        const embedding::MeshQuery empty_query(empty_mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        REQUIRE(!empty_query.BoundingBox().has_value());
        CHECK(!empty_query.IsWithinBoundingBox(PointType{}));
    }

    TEST_CASE("MeshQueryTestSuite - Rejects Invalid Ray Directions")
    {
        CHECK_THROWS_AS((void)Ray(PointType{}, Vector3d{}), queso::Exception);
        CHECK_THROWS_AS(
            (void)Ray(PointType{}, Vector3d{ std::numeric_limits<double>::infinity(), 0.0, 0.0 }), queso::Exception
        );
    }

    TEST_CASE("MeshQueryTestSuite - Rejects Invalid Aabb Diagnostics And Debug Candidate Preconditions")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        const PointType lower{ 0.0, 0.0, 0.0 };
        const PointType upper{ 1.0, 1.0, 1.0 };
        const std::vector<IndexType> invalid_candidates{ mesh.NumOfTriangles() };
        if constexpr (!NOTDEBUG) {
            CHECK_THROWS_AS((void)query.IntersectsAabb(invalid_candidates, lower, upper), Exception);
            CHECK_THROWS_AS((void)query.GetIntersectedTriangleIds(invalid_candidates, lower, upper), Exception);
        }
        if (!QuESo_DIAGNOSTICS_ENABLED()) { return; }
        const PointType non_finite{ std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0 };
        CHECK_THROWS_AS((void)query.GetAabbCandidates(non_finite, upper), Exception);
        CHECK_THROWS_AS((void)query.GetIntersectedTriangleIds({}, upper, lower), Exception);
        const double snap_distance = MakeTolerance().SnapDistance();
        const PointType thin_upper{ snap_distance, snap_distance, snap_distance };
        CHECK_THROWS_AS(
            (void)query.IntersectsAabb({}, lower, thin_upper, AabbIntersectionPolicy::SnapEroded), Exception
        );
        const PointType inset_boundary{ 2.0 * snap_distance, 2.0 * snap_distance, 2.0 * snap_distance };
        CHECK_NOTHROW((void)query.IntersectsAabb({}, lower, inset_boundary, AabbIntersectionPolicy::SnapEroded));
    }

    TEST_CASE("MeshQueryTestSuite - Direction Magnitude Does Not Change Classification")
    {
        const TriangleMesh mesh = MeshUtilities::MakeMeshBox({ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
        const embedding::MeshQuery query(mesh.View(), embedding::MeshQueryMode::Closed, MakeTolerance());
        constexpr PointType origin{ 0.5, 0.5, 0.5 };
        const auto unit_scale = query.Classify(Ray(origin, Vector3d{ 1.0, 0.123, 0.312 }));
        const auto large_scale = query.Classify(Ray(origin, Vector3d{ 10.0, 1.23, 3.12 }));
        CHECK((unit_scale.is_inside) == (large_scale.is_inside));
        CHECK((unit_scale.is_conclusive) == (large_scale.is_conclusive));
    }

}  // TEST_SUITE

}  // namespace queso::testing

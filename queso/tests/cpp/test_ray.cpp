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

//// Project includes
#include "helpers/numerical_checks.hpp"
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/ray.h"

namespace queso::testing {
namespace {

    using embedding::Ray;
    using embedding::RayTriangleIntersectionStatus;

    [[nodiscard]] TriangleMesh MakeTriangle(double Scale = 1.0)
    {
        TriangleMesh mesh;
        const IndexType first = mesh.AddVertex({ 0.0, 0.0, 0.0 });
        const IndexType second = mesh.AddVertex({ Scale, 0.0, 0.0 });
        const IndexType third = mesh.AddVertex({ 0.0, Scale, 0.0 });
        mesh.AddTriangle({ first, second, third }, { 0.0, 0.0, 1.0 });
        return mesh;
    }

}  // namespace

TEST_SUITE("RayTestSuite")
{

    TEST_CASE("RayTestSuite - Triangle Intersection Reports Named Statuses")
    {
        const TriangleMesh mesh = MakeTriangle();
        const auto triangle = mesh.Triangle<WithoutNormals>(0);
        const GeometryTolerance tolerance =
            GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

        const Ray hit_ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
        const auto hit = hit_ray.IntersectTriangle(triangle, tolerance.ZeroLength());
        CHECK(hit.status == RayTriangleIntersectionStatus::Hit);
        CheckNear(hit.distance, 1.0, 1e-12);

        const Ray miss_ray(PointType{ 2.0, 2.0, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
        CHECK(
            miss_ray.IntersectTriangle(triangle, tolerance.ZeroLength()).status == RayTriangleIntersectionStatus::Miss
        );

        const Ray parallel_ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.0, 0.0 });
        CHECK(
            parallel_ray.IntersectTriangle(triangle, tolerance.ZeroLength()).status
            == RayTriangleIntersectionStatus::Parallel
        );

        TriangleMesh degenerate_mesh;
        const IndexType first = degenerate_mesh.AddVertex({ 0.0, 0.0, 0.0 });
        const IndexType second = degenerate_mesh.AddVertex({ 1.0, 0.0, 0.0 });
        const IndexType third = degenerate_mesh.AddVertex({ 2.0, 0.0, 0.0 });
        degenerate_mesh.AddTriangle({ first, second, third }, { 0.0, 0.0, 0.0 });
        CHECK(
            hit_ray.IntersectTriangle(degenerate_mesh.Triangle<WithoutNormals>(0), tolerance.ZeroLength()).status
            == RayTriangleIntersectionStatus::Degenerate
        );
    }

    TEST_CASE("RayTestSuite - Parallel Classification Is Scale Independent")
    {
        const TriangleMesh unit_mesh = MakeTriangle();
        const TriangleMesh small_mesh = MakeTriangle(1e-6);
        const Ray parallel_ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.0, 0.0 });
        const Ray normal_ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
        for (const auto triangle : { unit_mesh.Triangle<WithoutNormals>(0), small_mesh.Triangle<WithoutNormals>(0) }) {
            CHECK(parallel_ray.IsParallel(triangle));
            CHECK(!normal_ray.IsParallel(triangle));
        }
    }

}  // TEST_SUITE

}  // namespace queso::testing

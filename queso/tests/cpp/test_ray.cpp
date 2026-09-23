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

//// Project includes
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/containers/triangle_mesh.hpp"
#include "queso/embedding/ray.h"
#include "queso/includes/checks.hpp"

namespace queso::Testing {
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

BOOST_AUTO_TEST_SUITE(RayTestSuite)

BOOST_AUTO_TEST_CASE(TriangleIntersectionReportsNamedStatuses)
{
    const TriangleMesh mesh = MakeTriangle();
    const auto triangle = mesh.Triangle<WithoutNormals>(0);
    const GeometryTolerance tolerance = GeometryTolerance::FromScale({ .length_scale = 1.0, .coordinate_scale = 1.0 });

    const Ray hit_ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
    const auto hit = hit_ray.IntersectTriangle(triangle, tolerance.ZeroLength());
    QuESo_CHECK(hit.status == RayTriangleIntersectionStatus::Hit);
    QuESo_CHECK_NEAR(hit.distance, 1.0, 1e-12);

    const Ray miss_ray(PointType{ 2.0, 2.0, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
    QuESo_CHECK(
        miss_ray.IntersectTriangle(triangle, tolerance.ZeroLength()).status == RayTriangleIntersectionStatus::Miss
    );

    const Ray parallel_ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.0, 0.0 });
    QuESo_CHECK(
        parallel_ray.IntersectTriangle(triangle, tolerance.ZeroLength()).status
        == RayTriangleIntersectionStatus::Parallel
    );

    TriangleMesh degenerate_mesh;
    const IndexType first = degenerate_mesh.AddVertex({ 0.0, 0.0, 0.0 });
    const IndexType second = degenerate_mesh.AddVertex({ 1.0, 0.0, 0.0 });
    const IndexType third = degenerate_mesh.AddVertex({ 2.0, 0.0, 0.0 });
    degenerate_mesh.AddTriangle({ first, second, third }, { 0.0, 0.0, 0.0 });
    QuESo_CHECK(
        hit_ray.IntersectTriangle(degenerate_mesh.Triangle<WithoutNormals>(0), tolerance.ZeroLength()).status
        == RayTriangleIntersectionStatus::Degenerate
    );
}

BOOST_AUTO_TEST_CASE(ParallelClassificationIsScaleIndependent)
{
    const TriangleMesh unit_mesh = MakeTriangle();
    const TriangleMesh small_mesh = MakeTriangle(1e-6);
    const Ray parallel_ray(PointType{ 0.2, 0.2, 0.0 }, Vector3d{ 1.0, 0.0, 0.0 });
    const Ray normal_ray(PointType{ 0.2, 0.2, -1.0 }, Vector3d{ 0.0, 0.0, 1.0 });
    for (const auto triangle : { unit_mesh.Triangle<WithoutNormals>(0), small_mesh.Triangle<WithoutNormals>(0) }) {
        QuESo_CHECK(parallel_ray.IsParallel(triangle));
        QuESo_CHECK(!normal_ray.IsParallel(triangle));
    }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace queso::Testing

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
#include <cmath>
#include <utility>

//// Project includes
#include "queso/embedding/local_surface_classifier.h"
#include "queso/embedding/ray.h"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso::embedding::detail {
namespace {

    [[nodiscard]] std::pair<bool, bool> ClassifyRay(const Ray& rRay, const TriangleMeshView& rMesh, double ZeroLength)
    {
        double min_distance = MAXD;
        bool is_inside = false;

        for (const auto& r_triangle : rMesh.Triangles<WithoutNormals>()) {
            const auto intersection = rRay.IntersectTriangle(r_triangle, ZeroLength);
            if (intersection.status != RayTriangleIntersectionStatus::Hit) { continue; }

            const double sum_u_v = intersection.u + intersection.v;
            if (intersection.distance < ZeroLength) { return { false, true }; }
            if (intersection.u < ray_predicates::BarycentricBoundaryTolerance
                || intersection.v < ray_predicates::BarycentricBoundaryTolerance
                || sum_u_v > 1.0 - ray_predicates::BarycentricBoundaryTolerance) {
                return { false, false };
            }
            if (intersection.distance < min_distance) {
                is_inside = intersection.is_back_facing;
                min_distance = intersection.distance;
            }
        }
        return { is_inside, true };
    }

}  // namespace

LocalSurfaceClassification
    ClassifyOnBoundedSide(const PointType& rPoint, const TriangleMeshView& rBoundarySection, double ZeroLength)
{
    const IndexType num_triangles = rBoundarySection.NumOfTriangles();
    if (num_triangles == 0) { return LocalSurfaceClassification::Inconclusive; }

    IndexType current_id = 0;
    IndexType success_count = 0;
    int inside_count = 0;
    while (success_count < 10 && current_id < num_triangles) {
        const auto triangle = rBoundarySection.Triangle<WithoutNormals>(current_id++);
        Vector3d direction = TriangleUtilities::Center(triangle) - rPoint;
        if (Math::SquaredNorm(direction) <= ZeroLength * ZeroLength) { continue; }

        const Ray ray(rPoint, direction);
        if (ray.IsParallel(triangle)) { continue; }

        const auto [is_inside, success] = ClassifyRay(ray, rBoundarySection, ZeroLength);
        if (success) {
            ++success_count;
            inside_count += is_inside ? 1 : -1;
        }
    }
    if (success_count == 0 || inside_count == 0) { return LocalSurfaceClassification::Inconclusive; }
    return inside_count > 0 ? LocalSurfaceClassification::Inside : LocalSurfaceClassification::Outside;
}

}  // namespace queso::embedding::detail

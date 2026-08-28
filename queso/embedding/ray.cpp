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

//// Own include
#include "queso/embedding/ray.h"

//// STL includes
#include <algorithm>
#include <array>
#include <cmath>

//// Project includes
#include "queso/includes/numerical_guards.hpp"
#include "queso/utilities/math_utilities.hpp"

namespace queso::embedding {

namespace {

    constexpr double RayAngularTolerance = 1e-13;

}  // namespace

Ray::Ray(PointView Origin, PointView Direction) : mOrigin{ Origin[0], Origin[1], Origin[2] }
{
    const double norm = std::hypot(Direction[0], Direction[1], Direction[2]);
    QuESo_ERROR_IF(
        !numerical_guards::IsSafeDivisor(norm) || !std::isfinite(Origin[0]) || !std::isfinite(Origin[1])
        || !std::isfinite(Origin[2])
    ) << "Ray origin and direction must be finite, and direction must be nonzero.\n";
    mDirection = PointType{ Direction[0] / norm, Direction[1] / norm, Direction[2] / norm };
    for (IndexType axis = 0; axis < 3; ++axis) { mInverseDirection[axis] = 1.0 / mDirection[axis]; }

    mDirectionIsStrictlyPositive = mDirection[0] > 0.0 && mDirection[1] > 0.0 && mDirection[2] > 0.0;
    if (!mDirectionIsStrictlyPositive) {
        for (IndexType axis = 0; axis < 3; ++axis) { mDirectionSigns[axis] = mInverseDirection[axis] < 0.0; }
    }
}

bool Ray::IntersectsAabb(PointView rLowerBound, PointView rUpperBound) const noexcept
{
    if (mDirectionIsStrictlyPositive) {
        // More efficient, but expects ray direction to be in positive direction: x>0, y>0, z>0.
        return IntersectsPositiveAabb(rLowerBound, rUpperBound);
    } else {
        /// Works with all Ray directions.
        return IntersectsGeneralAabb(rLowerBound, rUpperBound);
    }
}

bool Ray::IntersectsPositiveAabb(PointView rLowerBound, PointView rUpperBound) const noexcept
{
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    double lower_0 = rLowerBound[0];
    double lower_1 = rLowerBound[1];

    double upper_0 = rUpperBound[0];
    double upper_1 = rUpperBound[1];

    double origin_0 = mOrigin[0];
    double origin_1 = mOrigin[1];

    double inv_direction_0 = mInverseDirection[0];
    double inv_direction_1 = mInverseDirection[1];

    double lower_2 = rLowerBound[2];
    double upper_2 = rUpperBound[2];
    double origin_2 = mOrigin[2];

    // Check if origin of lies inside aabb.
    if (origin_0 >= lower_0 && origin_1 >= lower_1 && origin_2 >= lower_2 && origin_0 <= upper_0 && origin_1 <= upper_1
        && origin_2 <= upper_2) {
        return true;
    }

    tmin = (lower_0 - origin_0) * inv_direction_0;
    tymax = (upper_1 - origin_1) * inv_direction_1;
    if (tmin > tymax) { return false; }

    tmax = (upper_0 - origin_0) * inv_direction_0;
    tymin = (lower_1 - origin_1) * inv_direction_1;


    if (tymin > tmax) return false;

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    double inv_direction_2 = mInverseDirection[2];
    tzmin = (lower_2 - origin_2) * inv_direction_2;
    if ((tzmin > tmax)) return false;

    tzmax = (upper_2 - origin_2) * inv_direction_2;
    if (tmin > tzmax) return false;

    if (tzmin > tmin) tmin = tzmin;

    if (tmin < 0) { return false; }

    return true;
}

bool Ray::IntersectsGeneralAabb(PointView rLowerBound, PointView rUpperBound) const noexcept
{
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    std::array<Vector3d, 2> bounds = { Vector3d{ rLowerBound[0], rLowerBound[1], rLowerBound[2] },
                                       Vector3d{ rUpperBound[0], rUpperBound[1], rUpperBound[2] } };
    double lower_0 = bounds[mDirectionSigns[0]][0];
    double lower_1 = bounds[mDirectionSigns[1]][1];

    double upper_0 = bounds[1 - mDirectionSigns[0]][0];
    double upper_1 = bounds[1 - mDirectionSigns[1]][1];

    double origin_0 = mOrigin[0];
    double origin_1 = mOrigin[1];

    double inv_direction_0 = mInverseDirection[0];
    double inv_direction_1 = mInverseDirection[1];

    double lower_2 = bounds[mDirectionSigns[2]][2];
    double upper_2 = bounds[1 - mDirectionSigns[2]][2];
    double origin_2 = mOrigin[2];

    // Check if origin of lies inside aabb.
    if (origin_0 >= rLowerBound[0] && origin_1 >= rLowerBound[1] && origin_2 >= rLowerBound[2]
        && origin_0 <= rUpperBound[0] && origin_1 <= rUpperBound[1] && origin_2 <= rUpperBound[2]) {
        return true;
    }

    tmin = (lower_0 - origin_0) * inv_direction_0;
    tymax = (upper_1 - origin_1) * inv_direction_1;
    if (tmin > tymax) { return false; }

    tmax = (upper_0 - origin_0) * inv_direction_0;
    tymin = (lower_1 - origin_1) * inv_direction_1;


    if (tymin > tmax) return false;

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    double inv_direction_2 = mInverseDirection[2];
    tzmin = (lower_2 - origin_2) * inv_direction_2;
    if ((tzmin > tmax)) return false;

    tzmax = (upper_2 - origin_2) * inv_direction_2;
    if (tmin > tzmax) return false;

    if (tzmin > tmin) tmin = tzmin;

    if (tmin < 0) { return false; }

    return true;
}

RayTriangleIntersection
    Ray::IntersectTriangle(TriangleProxy<WithoutNormals> Triangle, double MinimumRayDistance) const noexcept
{
    const PointView v0 = Triangle.P1;
    const PointView v1 = Triangle.P2;
    const PointView v2 = Triangle.P3;
    const Vector3d v0v1 = v1 - v0;
    const Vector3d v0v2 = v2 - v0;
    const Vector3d normal = Math::Cross(v0v1, v0v2);
    const double squared_normal_norm = Math::SquaredNorm(normal);
    if (squared_normal_norm == 0.0) { return { .status = RayTriangleIntersectionStatus::Degenerate }; }

    // Cross product: mDirection x v0v2
    const Vector3d pvec = Math::Cross(mDirection, v0v2);

    // Dot product: v0v1 * pvec
    const double det = Math::Dot(v0v1, pvec);
    if (det * det <= RayAngularTolerance * RayAngularTolerance * squared_normal_norm) {
        return { .status = RayTriangleIntersectionStatus::Parallel };
    }

    RayTriangleIntersection result{ .is_back_facing = det < 0.0 };

    // Get inverse of determinant.
    const double invDet = 1.0 / det;

    // Substraction: mOrigin - v0
    const Vector3d tvec = mOrigin - v0;

    // Dot product x invDet: (tvec * pvec) * invDet
    result.u = Math::Dot(tvec, pvec) * invDet;

    if (result.u < -detail::ray_predicates::BarycentricBoundaryTolerance
        || result.u > 1.0 + detail::ray_predicates::BarycentricBoundaryTolerance) {
        return result;
    }

    // Cross product: tvec x v0v1
    const Vector3d qvec = Math::Cross(tvec, v0v1);

    // Dot product x invDet: (mDirection * qvec) * invDet
    result.v = Math::Dot(mDirection, qvec) * invDet;

    if (result.v < -detail::ray_predicates::BarycentricBoundaryTolerance
        || result.u + result.v > 1.0 + detail::ray_predicates::BarycentricBoundaryTolerance) {
        return result;
    }

    // Dot product x invDet: (v0v2 * qvec) * invDet
    result.distance = Math::Dot(v0v2, qvec) * invDet;

    // Return false if ray intersects in negative direction.
    if (result.distance < -MinimumRayDistance) { return result; }

    result.status = RayTriangleIntersectionStatus::Hit;
    return result;
}

bool Ray::IsParallel(TriangleProxy<WithoutNormals> Triangle) const noexcept(NOTDEBUG)
{
    const Vector3d normal = Math::Cross(Triangle.P2 - Triangle.P1, Triangle.P3 - Triangle.P1);
    const double squared_normal_norm = Math::SquaredNorm(normal);
    if (squared_normal_norm == 0.0) { return true; }
    const double normal_component = Math::Dot(mDirection, normal);
    return normal_component * normal_component <= RayAngularTolerance * RayAngularTolerance * squared_normal_norm;
}

}  // namespace queso::embedding

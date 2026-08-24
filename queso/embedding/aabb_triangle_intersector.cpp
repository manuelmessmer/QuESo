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
#include "queso/embedding/aabb_triangle_intersector.h"

//// STL includes
#include <algorithm>
#include <cmath>

//// Project includes
#include "queso/utilities/math_utilities.hpp"

namespace queso::embedding::detail {

bool AabbTriangleIntersector::IntersectsTriangle(TriangleProxy<WithoutNormals> Triangle, double Inset) const noexcept
{
    const PointView v0 = Triangle.P1;
    const PointView v1 = Triangle.P2;
    const PointView v2 = Triangle.P3;

    // Erode the box uniformly before applying the separating-axis tests.
    const Vector3d extent{ (upperBound[0] - lowerBound[0]) / 2.0 - Inset,
                           (upperBound[1] - lowerBound[1]) / 2.0 - Inset,
                           (upperBound[2] - lowerBound[2]) / 2.0 - Inset };

    const Vector3d v0_orig{ v0[0] - centre[0], v0[1] - centre[1], v0[2] - centre[2] };
    const Vector3d v1_orig{ v1[0] - centre[0], v1[1] - centre[1], v1[2] - centre[2] };
    const Vector3d v2_orig{ v2[0] - centre[0], v2[1] - centre[1], v2[2] - centre[2] };

    // Compute the edge vectors of the triangle  (ABC). Line between vertices.
    const Vector3d f0 = v1 - v0;
    const Vector3d f1 = v2 - v1;
    const Vector3d f2 = v0 - v2;

    // Axis-aligned box basis vectors.
    const Vector3d u0{ 1.0, 0.0, 0.0 };
    const Vector3d u1{ 0.0, 1.0, 0.0 };
    const Vector3d u2{ 0.0, 0.0, 1.0 };

    // Test 13 separating-axis candidates: nine edge cross-products, three box axes, and the triangle normal.

    // First test (u0, u1, u2) vs. (f0, f1, f2). 9 tests in total.
    // u0 vs f0.
    Vector3d axis_u0_f0{ u0[1] * f0[2] - u0[2] * f0[1], u0[2] * f0[0] - u0[0] * f0[2], u0[0] * f0[1] - u0[1] * f0[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f0)) { return false; }

    // u0 vs f1.
    Vector3d axis_u0_f1{ u0[1] * f1[2] - u0[2] * f1[1], u0[2] * f1[0] - u0[0] * f1[2], u0[0] * f1[1] - u0[1] * f1[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f1)) { return false; }

    // u0 vs f2.
    Vector3d axis_u0_f2{ u0[1] * f2[2] - u0[2] * f2[1], u0[2] * f2[0] - u0[0] * f2[2], u0[0] * f2[1] - u0[1] * f2[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f2)) { return false; }

    // u1 vs f0.
    Vector3d axis_u1_f0{ u1[1] * f0[2] - u1[2] * f0[1], u1[2] * f0[0] - u1[0] * f0[2], u1[0] * f0[1] - u1[1] * f0[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f0)) { return false; }

    // u1 vs f1.
    Vector3d axis_u1_f1{ u1[1] * f1[2] - u1[2] * f1[1], u1[2] * f1[0] - u1[0] * f1[2], u1[0] * f1[1] - u1[1] * f1[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f1)) { return false; }

    // u1 vs f2.
    Vector3d axis_u1_f2{ u1[1] * f2[2] - u1[2] * f2[1], u1[2] * f2[0] - u1[0] * f2[2], u1[0] * f2[1] - u1[1] * f2[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f2)) { return false; }

    // u2 vs f0.
    Vector3d axis_u2_f0{ u2[1] * f0[2] - u2[2] * f0[1], u2[2] * f0[0] - u2[0] * f0[2], u2[0] * f0[1] - u2[1] * f0[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f0)) { return false; }

    // u2 vs f1.
    Vector3d axis_u2_f1{ u2[1] * f1[2] - u2[2] * f1[1], u2[2] * f1[0] - u2[0] * f1[2], u2[0] * f1[1] - u2[1] * f1[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f1)) { return false; }

    // u2 vs f2.
    Vector3d axis_u2_f2{ u2[1] * f2[2] - u2[2] * f2[1], u2[2] * f2[0] - u2[0] * f2[2], u2[0] * f2[1] - u2[1] * f2[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f2)) { return false; }

    // Test face normals of aabb. 3 Tests.
    // axis1: (1, 0, 0)
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u0)) { return false; }
    // axis2: (0, 1, 0)
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u1)) { return false; }

    // axis3 (0, 0, 1)
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u2)) { return false; }

    // Test face normal of triangle
    Vector3d triangle_normal{ f0[1] * f1[2] - f0[2] * f1[1],
                              f0[2] * f1[0] - f0[0] * f1[2],
                              f0[0] * f1[1] - f0[1] * f1[0] };
    if (!CheckAxis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, triangle_normal)) { return false; }

    // No separating axis was found.
    return true;
}


bool AabbTriangleIntersector::CheckAxis(
    PointView rAxisX,
    PointView rAxisY,
    PointView rAxisZ,
    PointView rFirst,
    PointView rSecond,
    PointView rThird,
    PointView rExtent,
    PointView rTestAxis
) const noexcept
{
    const double first_projection = Math::Dot(rFirst, rTestAxis);
    const double second_projection = Math::Dot(rSecond, rTestAxis);
    const double third_projection = Math::Dot(rThird, rTestAxis);
    const double radius = rExtent[0] * std::abs(Math::Dot(rAxisX, rTestAxis))
                          + rExtent[1] * std::abs(Math::Dot(rAxisY, rTestAxis))
                          + rExtent[2] * std::abs(Math::Dot(rAxisZ, rTestAxis));
    if (std::max(
            -std::max({ first_projection, second_projection, third_projection }),
            std::min({ first_projection, second_projection, third_projection })
        )
        > radius) {
        return false;
    }

    return true;
}

}  // namespace queso::embedding::detail

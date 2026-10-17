/*
  ____        ______  _____
 / __ \      |  ____|/ ____|
| |  | |_   _| |__  | (___   ___
| |  | | | | |  __|  \___ \ / _ \
| |__| | |_| | |____ ____) | (_) |
 \___\_\\__,_|______|_____/ \___/
        Quadrature for Embedded Solids

 License:    BSD 4-Clause License
             See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE

 Authors:    Manuel Messmer
*/

#pragma once

//// STL includes
#include <algorithm>
#include <concepts>
#include <functional>
#include <limits>
#include <numbers>

//// Project includes
#include "queso/containers/integration_point_concepts.hpp"
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/define.hpp"
#include "queso/includes/numerical_guards.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_gauss_legendre_integration_points.hpp"

namespace queso::triangle_utilities {

/// @brief Returns the physical triangle area in squared coordinate units.
/// @tparam TMode Triangle proxy mode, with or without stored normals.
/// @param rTriangle Triangle proxy.
/// @return Non-negative area. Zero for a degenerate triangle.
template<class TMode>
[[nodiscard]] inline double Area(const TriangleProxy<TMode>& rTriangle) noexcept
{
    const Vector3d edge_a = rTriangle.P2 - rTriangle.P1;
    const Vector3d edge_b = rTriangle.P3 - rTriangle.P1;

    return 0.5 * Math::Norm(Math::Cross(edge_a, edge_b));
}

/// @brief Computes the unit normal from vertex winding for a triangle without stored normals.
/// @param rTriangle Triangle proxy without a normal.
/// @return Unit normal, or zero if its magnitude is not a safe divisor.
[[nodiscard]] inline Vector3d Normal(const TriangleProxy<WithoutNormals>& rTriangle) noexcept
{
    const Vector3d edge_a = rTriangle.P2 - rTriangle.P1;
    const Vector3d edge_b = rTriangle.P3 - rTriangle.P2;
    const Vector3d edge_c = rTriangle.P1 - rTriangle.P3;

    const double length_a = Math::Norm(edge_a);
    const double length_b = Math::Norm(edge_b);
    const double length_c = Math::Norm(edge_c);

    PointType result{};
    if (length_a >= length_c && length_b >= length_c) {
        result = Math::Cross(edge_a, edge_b);
    } else if (length_a >= length_b && length_c >= length_b) {
        result = Math::Cross(edge_c, edge_a);
    } else {
        result = Math::Cross(edge_b, edge_c);
    }

    const double magnitude = Math::Norm(result);
    if (numerical_guards::IsSafeDivisor(magnitude)) {
        result *= 1.0 / magnitude;
    } else {
        result = { 0.0, 0.0, 0.0 };
    }
    return result;
}

/// @brief Returns the longest-edge aspect ratio, normalized to one for an equilateral triangle.
/// @tparam TMode Triangle proxy mode, with or without stored normals.
/// @param rTriangle Triangle proxy.
/// @return Dimensionless aspect ratio, or infinity if the area is not a safe divisor.
template<class TMode>
[[nodiscard]] inline double AspectRatio(const TriangleProxy<TMode>& rTriangle) noexcept
{
    const double triangle_area = Area(rTriangle);

    const double a = Math::Norm(rTriangle.P2 - rTriangle.P1);
    const double b = Math::Norm(rTriangle.P3 - rTriangle.P2);
    const double c = Math::Norm(rTriangle.P1 - rTriangle.P3);

    const double max_edge = std::max({ a, b, c });

    if (!numerical_guards::IsSafeDivisor(triangle_area)) { return std::numeric_limits<double>::infinity(); }

    return max_edge * (a + b + c) / (4.0 * std::numbers::sqrt3_v<double> * triangle_area);
}

/// @brief Returns the triangle centroid in the input coordinate system.
/// @tparam TMode Triangle proxy mode, with or without stored normals.
/// @param rTriangle Triangle proxy.
/// @return Arithmetic mean of the three vertices; supports constant evaluation.
template<class TMode>
[[nodiscard]] constexpr Vector3d Center(const TriangleProxy<TMode>& rTriangle) noexcept
{
    const auto& p1 = rTriangle.P1;
    const auto& p2 = rTriangle.P2;
    const auto& p3 = rTriangle.P3;

    return { 1.0 / 3.0 * (p1[0] + p2[0] + p3[0]),
             1.0 / 3.0 * (p1[1] + p2[1] + p3[1]),
             1.0 / 3.0 * (p1[2] + p2[2] + p3[2]) };
}

/// @brief Visits physical boundary integration points without allocating an intermediate collection.
/// @tparam TBoundaryIntegrationPointType Boundary point type constructible from coordinates, weight, and normal.
/// @tparam TMode Triangle proxy mode, with or without stored normals.
/// @tparam TCallback Callable invoked as an lvalue with a temporary boundary integration point.
/// @param rTriangle Triangle proxy.
/// @param Order Maximum polynomial degree integrated exactly.
/// @param rCallback Invoked once per point, synchronously. Never retained.
/// @throws Exception If Order is unsupported, before invoking the callback. Callback exceptions propagate unchanged.
template<concepts::BoundaryIntegrationPoint TBoundaryIntegrationPointType, class TMode, class TCallback>
    requires std::constructible_from<TBoundaryIntegrationPointType, double, double, double, double, const Vector3d&>
             && std::invocable<TCallback&, TBoundaryIntegrationPointType>
void VisitBoundaryIps(const TriangleProxy<TMode>& rTriangle, TriangleQuadratureOrder Order, TCallback&& rCallback)
{
    const auto rule = detail::GetRule(Order);
    const double triangle_area = Area(rTriangle);
    const Vector3d triangle_normal = [&] {
        if constexpr (std::same_as<TMode, WithNormals>) {
            return Vector3d{ rTriangle.Normal[0], rTriangle.Normal[1], rTriangle.Normal[2] };
        } else {
            return Normal(rTriangle);
        }
    }();

    for (SizeType i = 0; i < rule.points.size(); ++i) {
        const auto& shape_values = rule.shape_functions[i];
        PointType point{};
        for (IndexType d = 0; d < 3; ++d) {
            point[d] = shape_values[0] * rTriangle.P1[d] + shape_values[1] * rTriangle.P2[d]
                       + shape_values[2] * rTriangle.P3[d];
        }
        const double weight = 2.0 * rule.points[i].Weight() * triangle_area;
        std::invoke(rCallback, TBoundaryIntegrationPointType{ point[0], point[1], point[2], weight, triangle_normal });
    }
}

}  // namespace queso::triangle_utilities

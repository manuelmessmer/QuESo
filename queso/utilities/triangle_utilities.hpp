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

//// Project includes
#include <numbers>

//// Project includes
#include "queso/containers/integration_point.hpp"
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/define.hpp"
#include "queso/includes/numerical_guards.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_gauss_legendre_integration_points.hpp"

namespace queso {
namespace TriangleUtilities {

    /// @brief Returns Area of triangle.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return double
    template<class Mode>
    inline double Area(const TriangleProxy<Mode>& rTriangle)
    {
        const Vector3d edge_a = rTriangle.P2 - rTriangle.P1;
        const Vector3d edge_b = rTriangle.P3 - rTriangle.P1;

        return 0.5 * Math::Norm(Math::Cross(edge_a, edge_b));
    }

    /// @brief Returns normal computed via vertices.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return Vector3d
    template<class Mode>
    inline Vector3d Normal(const TriangleProxy<Mode>& rTriangle)
    {
        const Vector3d edge_a = rTriangle.P2 - rTriangle.P1;
        const Vector3d edge_b = rTriangle.P3 - rTriangle.P2;
        const Vector3d edge_c = rTriangle.P1 - rTriangle.P3;

        const double length_a = Math::Norm(edge_a);
        const double length_b = Math::Norm(edge_b);
        const double length_c = Math::Norm(edge_c);

        PointType normal{};
        if (length_a >= length_c && length_b >= length_c) {
            normal = Math::Cross(edge_a, edge_b);
        } else if (length_a >= length_b && length_c >= length_b) {
            normal = Math::Cross(edge_c, edge_a);
        } else {
            normal = Math::Cross(edge_b, edge_c);
        }

        const double norm = Math::Norm(normal);
        if (numerical_guards::IsSafeDivisor(norm)) {
            normal *= 1.0 / norm;
        } else {
            normal = { 0.0, 0.0, 0.0 };
        }
        return normal;
    }

    /// @brief Returns AspectRatio of triangle.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return double
    template<class Mode>
    inline double AspectRatio(const TriangleProxy<Mode>& rTriangle)
    {
        const auto area = Area(rTriangle);

        const double a = Math::Norm(rTriangle.P2 - rTriangle.P1);  // length a
        const double b = Math::Norm(rTriangle.P3 - rTriangle.P2);  // length b
        const double c = Math::Norm(rTriangle.P1 - rTriangle.P3);  // length c

        const double max_edge = std::max({ a, b, c });

        if (!numerical_guards::IsSafeDivisor(area)) { return std::numeric_limits<double>::infinity(); }

        return max_edge * (a + b + c) / (4.0 * std::numbers::sqrt3_v<double> * area);
    }

    /// @brief Center of triangles in global coordinates.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return Vector3d.
    template<class Mode>
    inline Vector3d Center(const TriangleProxy<Mode>& rTriangle)
    {
        const auto& p1 = rTriangle.P1;
        const auto& p2 = rTriangle.P2;
        const auto& p3 = rTriangle.P3;

        return { 1.0 / 3.0 * (p1[0] + p2[0] + p3[0]),
                 1.0 / 3.0 * (p1[1] + p2[1] + p3[1]),
                 1.0 / 3.0 * (p1[2] + p2[2] + p3[2]) };
    }

    namespace detail {
        using IpVectorType = std::vector<IntegrationPoint>;
        inline const IpVectorType& GetIntegrationPoints(IndexType method)
        {
            switch (method) {
            case 0:
                return TriangleGaussLegendrePoints1::IntegrationPoints();
            case 1:
                return TriangleGaussLegendrePoints2::IntegrationPoints();
            case 2:
                return TriangleGaussLegendrePoints3::IntegrationPoints();
            case 3:
                return TriangleGaussLegendrePoints4::IntegrationPoints();
            case 4:
                return TriangleGaussLegendrePoints4::IntegrationPoints();
            }
            QuESo_ERROR << "Wrong Index of Shape Function.\n";
        }

        inline double ShapeFunctionValue(IndexType ShapeFunctionIndex, const Vector3d& rPoint)
        {
            switch (ShapeFunctionIndex) {
            case 0:
                return (1.0 - rPoint[0] - rPoint[1]);
            case 1:
                return (rPoint[0]);
            case 2:
                return (rPoint[1]);
            }
            QuESo_ERROR << "Wrong Index of Shape Function.\n";
        }
    }  // namespace detail

    /// @brief Get boundary integration points in global space.
    /// @tparam TBoundaryIntegrationPointType
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @param Method integration method.
    /// @return Boundary integration points.
    /// TODO: refactor to VisitBoundaryIps with TCallback to direcly push_back to a vector (with or without
    /// transformation).
    /// TriangleUtilities::VisitIPsGlobal<TBoundaryIntegrationPointType>(rTriangle, method, [&](auto&& rIp) {
    ///     if constexpr (TSpace == CoordinateSpace::global) {
    ///         boundary_ips.push_back(std::forward<decltype(rIp)>(rIp));
    ///     } else {
    ///         boundary_ips.push_back(rCellMapper.ToParametric(rIp));
    ///     }
    /// });
    template<typename TBoundaryIntegrationPointType, class Mode>
    std::vector<TBoundaryIntegrationPointType> GetIPsGlobal(const TriangleProxy<Mode>& rTriangle, IndexType Method)
    {
        const auto& s_integration_points = detail::GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        auto global_integration_points = std::vector<TBoundaryIntegrationPointType>();

        const auto& p1 = rTriangle.P1;
        const auto& p2 = rTriangle.P2;
        const auto& p3 = rTriangle.P3;

        for (IndexType i = 0; i < point_numbers; ++i) {
            const double x = detail::ShapeFunctionValue(0, s_integration_points[i].Point()) * p1[0]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].Point()) * p2[0]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].Point()) * p3[0];

            const double y = detail::ShapeFunctionValue(0, s_integration_points[i].Point()) * p1[1]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].Point()) * p2[1]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].Point()) * p3[1];

            const double z = detail::ShapeFunctionValue(0, s_integration_points[i].Point()) * p1[2]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].Point()) * p2[2]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].Point()) * p3[2];

            // Normalize weights to 1 by multiplying by 2.
            const double weight = 2.0 * s_integration_points[i].Weight() * Area(rTriangle);
            global_integration_points.emplace_back(x, y, z, weight, Normal(rTriangle));
        }

        return global_integration_points;
    }

}  // namespace TriangleUtilities
}  // namespace queso

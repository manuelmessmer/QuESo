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
#include "queso/containers/integration_point.hpp"
#include "queso/containers/triangle_mesh_concepts.hpp"
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_gauss_legendre_integration_points.hpp"

namespace queso {
namespace TriangleUtilities {

    /// @brief Returns Area of triangle.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return double
    template<class Mode>
    inline double Area(const TriangleProxy<Mode> &rTriangle)
    {
        const Vector3d A = rTriangle.P2 - rTriangle.P1;
        const Vector3d B = rTriangle.P3 - rTriangle.P1;

        return 0.5 * Math::Norm(Math::Cross(A, B));
    }

    /// @brief Returns normal computed via vertices.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return Vector3d
    template<class Mode>
    inline Vector3d Normal(const TriangleProxy<Mode> &rTriangle)
    {
        const Vector3d A = rTriangle.P2 - rTriangle.P1;
        const Vector3d B = rTriangle.P3 - rTriangle.P2;
        const Vector3d C = rTriangle.P1 - rTriangle.P3;

        const double lenght_A = Math::Norm(A);
        const double lenght_B = Math::Norm(B);
        const double lenght_C = Math::Norm(C);

        PointType normal{};
        if (lenght_A >= lenght_C - ZEROTOL && lenght_B >= lenght_C - ZEROTOL) {
            normal = Math::Cross(A, B);
        } else if (lenght_A >= lenght_B - ZEROTOL && lenght_C >= lenght_B - ZEROTOL) {
            normal = Math::Cross(C, A);
        } else {
            normal = Math::Cross(B, C);
        }

        const double norm = Math::Norm(normal);
        if (norm > ZEROTOL) {
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
    inline double AspectRatio(const TriangleProxy<Mode> &rTriangle)
    {
        const auto area = Area(rTriangle);

        const double a = Math::Norm(rTriangle.P2 - rTriangle.P1);// length a
        const double b = Math::Norm(rTriangle.P3 - rTriangle.P2);// length b
        const double c = Math::Norm(rTriangle.P1 - rTriangle.P3);// length c

        const double max_edge = std::max(std::max(a, b), c);

        const double square_root_3 = 1.73205080756887729;
        if (area > EPS4) {
            return max_edge * (a + b + c) / (4.0 * square_root_3 * area);
        } else {
            return 1e10;
        }
    }

    /// @brief Center of triangles in global coordinates.
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @return Vector3d.
    template<class Mode>
    inline Vector3d Center(const TriangleProxy<Mode> &rTriangle)
    {
        const auto &P1 = rTriangle.P1;
        const auto &P2 = rTriangle.P2;
        const auto &P3 = rTriangle.P3;

        return { 1.0 / 3.0 * (P1[0] + P2[0] + P3[0]),
            1.0 / 3.0 * (P1[1] + P2[1] + P3[1]),
            1.0 / 3.0 * (P1[2] + P2[2] + P3[2]) };
    }

    namespace detail {
        using IpVectorType = std::vector<IntegrationPoint>;
        inline const IpVectorType &GetIntegrationPoints(IndexType method)
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

        inline double ShapeFunctionValue(IndexType ShapeFunctionIndex, const Vector3d &rPoint)
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
    }// namespace detail

    /// @brief Get boundary integration points in global space.
    /// @tparam TBoundaryIntegrationPointType
    /// @tparam Mode Triangle mode.
    /// @param rTriangle Triangle proxy.
    /// @param Method integration method.
    /// @return Boundary integration points.
    template<typename TBoundaryIntegrationPointType, class Mode>
    std::vector<TBoundaryIntegrationPointType> GetIPsGlobal(const TriangleProxy<Mode> &rTriangle, IndexType Method)
    {
        const auto &s_integration_points = detail::GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        auto global_integration_points = std::vector<TBoundaryIntegrationPointType>();

        const auto &P1 = rTriangle.P1;
        const auto &P2 = rTriangle.P2;
        const auto &P3 = rTriangle.P3;

        for (IndexType i = 0; i < point_numbers; ++i) {
            const double x = detail::ShapeFunctionValue(0, s_integration_points[i].data()) * P1[0]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].data()) * P2[0]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].data()) * P3[0];

            const double y = detail::ShapeFunctionValue(0, s_integration_points[i].data()) * P1[1]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].data()) * P2[1]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].data()) * P3[1];

            const double z = detail::ShapeFunctionValue(0, s_integration_points[i].data()) * P1[2]
                             + detail::ShapeFunctionValue(1, s_integration_points[i].data()) * P2[2]
                             + detail::ShapeFunctionValue(2, s_integration_points[i].data()) * P3[2];

            // Normalize weights to 1 by multiplying by 2.
            const double weight = 2.0 * s_integration_points[i].Weight() * Area(rTriangle);
            global_integration_points.emplace_back(x, y, z, weight, Normal(rTriangle));
        }

        return global_integration_points;
    }

}// namespace TriangleUtilities
}// namespace queso

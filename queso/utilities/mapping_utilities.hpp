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

#pragma once

//// STL includes
#include <concepts>

//// Project includes
#include "queso/containers/element_core.hpp"
#include "queso/containers/integration_point_concepts.hpp"
#include "queso/includes/define.hpp"
#include "queso/utilities/math_utilities.hpp"

namespace queso::mapping {

namespace detail {

    /// @brief Computes axis-wise affine scaling between two bounding boxes.
    /// @param rFrom Source bounds.
    /// @param rTo Target bounds.
    /// @return Vector3d
    [[nodiscard]] constexpr Vector3d
        ComputeScale(const BoundingBoxType& rFrom, const BoundingBoxType& rTo) noexcept(NOTDEBUG)
    {
        QuESo_ASSERT(std::abs(rFrom.upper[0] - rFrom.lower[0]) > ZEROTOL, "Degenerate bounding box.");
        QuESo_ASSERT(std::abs(rFrom.upper[1] - rFrom.lower[1]) > ZEROTOL, "Degenerate bounding box.");
        QuESo_ASSERT(std::abs(rFrom.upper[2] - rFrom.lower[2]) > ZEROTOL, "Degenerate bounding box.");
        return { std::abs(rTo.upper[0] - rTo.lower[0]) / std::abs(rFrom.upper[0] - rFrom.lower[0]),
                 std::abs(rTo.upper[1] - rTo.lower[1]) / std::abs(rFrom.upper[1] - rFrom.lower[1]),
                 std::abs(rTo.upper[2] - rTo.lower[2]) / std::abs(rFrom.upper[2] - rFrom.lower[2]) };
    }

    /// @brief Applies the Piola-type area-vector scaling induced by an affine box map.
    /// @param rNormal Input surface normal.
    /// @param rScale Axis-wise scale from source to target space.
    /// @return Vector3d
    [[nodiscard]] constexpr Vector3d TransformAreaVector(Vector3dView rNormal, const Vector3d& rScale) noexcept
    {
        return { rNormal[0] * rScale[1] * rScale[2],
                 rNormal[1] * rScale[0] * rScale[2],
                 rNormal[2] * rScale[0] * rScale[1] };
    }

    /// @brief Transforms a surface normal and returns the corresponding area scale.
    /// @param rNormal Input surface normal.
    /// @param rScale Axis-wise scale from source to target space.
    /// @return std::pair<Vector3d, double>
    [[nodiscard]] inline std::pair<Vector3d, double> TransformAreaData(Vector3dView rNormal, const Vector3d& rScale)
    {
        auto transformed = TransformAreaVector(rNormal, rScale);
        const double norm = Math::Norm(transformed);
        if (norm > ZEROTOL) {
            transformed *= 1.0 / norm;
        } else {
            transformed = { 0.0, 0.0, 0.0 };
        }
        return { transformed, norm };
    }

}// namespace detail

/// @brief Maps a point from global to parametric space using element-local bounds.
/// @param rPoint Point in global coordinates.
/// @param rBounds Element-local bounds.
/// @return PointType
[[nodiscard]] constexpr PointType ToParametric(PointView rPoint, const ElementBounds& rBounds) noexcept(NOTDEBUG)
{
    const auto scale = detail::ComputeScale(rBounds.global, rBounds.parametric);
    return { (rPoint[0] - rBounds.global.lower[0]) * scale[0] + rBounds.parametric.lower[0],
             (rPoint[1] - rBounds.global.lower[1]) * scale[1] + rBounds.parametric.lower[1],
             (rPoint[2] - rBounds.global.lower[2]) * scale[2] + rBounds.parametric.lower[2] };
}

/// @brief Maps a point from parametric to global space using element-local bounds.
/// @param rPoint Point in parametric coordinates.
/// @param rBounds Element-local bounds.
/// @return PointType
[[nodiscard]] constexpr PointType ToGlobal(PointView rPoint, const ElementBounds& rBounds) noexcept(NOTDEBUG)
{
    const auto scale = detail::ComputeScale(rBounds.parametric, rBounds.global);
    return { (rPoint[0] - rBounds.parametric.lower[0]) * scale[0] + rBounds.global.lower[0],
             (rPoint[1] - rBounds.parametric.lower[1]) * scale[1] + rBounds.global.lower[1],
             (rPoint[2] - rBounds.parametric.lower[2]) * scale[2] + rBounds.global.lower[2] };
}

/// @brief Returns the determinant of the element-local affine map Jacobian.
/// @param rBounds Element-local bounds.
/// @return double
[[nodiscard]] constexpr double DetJ(const ElementBounds& rBounds) noexcept(NOTDEBUG)
{
    const auto scale = detail::ComputeScale(rBounds.parametric, rBounds.global);
    return scale[0] * scale[1] * scale[2];
}

/// @brief Maps a volume integration point from global to parametric space.
/// @tparam TIp Integration-point type constructible from `(PointType, double)`.
/// @param rIp Integration point in global coordinates.
/// @param rBounds Element-local bounds.
/// @return TIp
template<concepts::IntegrationPoint TIp>
    requires(!concepts::BoundaryIntegrationPoint<TIp>)
[[nodiscard]] constexpr TIp ToParametric(const TIp& rIp, const ElementBounds& rBounds)
{
    static_assert(
        std::constructible_from<TIp, PointType, double>, "Integration point mapping requires TIp(PointType, double)."
    );
    return TIp(ToParametric(rIp.Point(), rBounds), rIp.Weight() / DetJ(rBounds));
}

/// @brief Maps a volume integration point from parametric to global space.
/// @tparam TIp Integration-point type constructible from `(PointType, double)`.
/// @param rIp Integration point in parametric coordinates.
/// @param rBounds Element-local bounds.
/// @return TIp
template<concepts::IntegrationPoint TIp>
    requires(!concepts::BoundaryIntegrationPoint<TIp>)
[[nodiscard]] constexpr TIp ToGlobal(const TIp& rIp, const ElementBounds& rBounds)
{
    static_assert(
        std::constructible_from<TIp, PointType, double>, "Integration point mapping requires TIp(PointType, double)."
    );
    return TIp(ToGlobal(rIp.Point(), rBounds), rIp.Weight() * DetJ(rBounds));
}

/// @brief Maps a boundary integration point from global to parametric space.
/// @tparam TBoundaryIp Boundary integration-point type constructible from `(PointType, double, Vector3d)`.
/// @param rIp Boundary integration point in global coordinates.
/// @param rBounds Element-local bounds.
/// @return TBoundaryIp
template<concepts::BoundaryIntegrationPoint TBoundaryIp>
[[nodiscard]] inline TBoundaryIp ToParametric(const TBoundaryIp& rIp, const ElementBounds& rBounds)
{
    static_assert(
        std::constructible_from<TBoundaryIp, PointType, double, Vector3d>,
        "Boundary integration point mapping requires TBoundaryIp(PointType, double, Vector3d)."
    );
    const auto [normal, area_scale] =
        detail::TransformAreaData(rIp.Normal(), detail::ComputeScale(rBounds.global, rBounds.parametric));
    return TBoundaryIp(ToParametric(rIp.Point(), rBounds), rIp.Weight() * area_scale, normal);
}

/// @brief Maps a boundary integration point from parametric to global space.
/// @tparam TBoundaryIp Boundary integration-point type constructible from `(PointType, double, Vector3d)`.
/// @param rIp Boundary integration point in parametric coordinates.
/// @param rBounds Element-local bounds.
/// @return TBoundaryIp
template<concepts::BoundaryIntegrationPoint TBoundaryIp>
[[nodiscard]] inline TBoundaryIp ToGlobal(const TBoundaryIp& rIp, const ElementBounds& rBounds)
{
    static_assert(
        std::constructible_from<TBoundaryIp, PointType, double, Vector3d>,
        "Boundary integration point mapping requires TBoundaryIp(PointType, double, Vector3d)."
    );
    const auto [normal, area_scale] =
        detail::TransformAreaData(rIp.Normal(), detail::ComputeScale(rBounds.parametric, rBounds.global));
    return TBoundaryIp(ToGlobal(rIp.Point(), rBounds), rIp.Weight() * area_scale, normal);
}

}// namespace queso::mapping

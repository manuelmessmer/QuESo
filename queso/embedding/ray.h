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
#include <cstdint>

//// Project includes
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/define.hpp"

namespace queso::embedding {

namespace detail::ray_predicates {

    inline constexpr double BarycentricBoundaryTolerance = 1e-14;

}  // namespace detail::ray_predicates

/// @brief Outcome category of a positive ray-triangle query.
enum class RayTriangleIntersectionStatus : std::uint8_t {
    Miss,  ///< The positive ray does not intersect the triangle interior or tolerance boundary.
    Hit,  ///< The positive ray intersects the triangle; all result fields are valid.
    Parallel,  ///< The ray direction is parallel to the non-degenerate triangle plane.
    Degenerate  ///< The triangle has exactly zero area.
};

/// @brief Result of intersecting a normalized positive ray with one triangle.
struct RayTriangleIntersection
{
    RayTriangleIntersectionStatus status = RayTriangleIntersectionStatus::Miss;  ///< Outcome category.
    double distance{};  ///< Physical distance from ray origin; meaningful only for Hit.
    double u{};  ///< First barycentric coordinate; meaningful only for Hit.
    double v{};  ///< Second barycentric coordinate; meaningful only for Hit.
    bool is_back_facing{};  ///< Whether the hit triangle is back-facing; meaningful only for Hit.
};

/// @brief Normalized ray supporting conservative box queries and exact triangle intersection.
class Ray
{
public:
    /// @brief Constructs a ray and normalizes its direction.
    /// @details Throws when the origin or direction is non-finite, or when the direction has zero length.
    /// @param Origin Finite ray origin.
    /// @param Direction Finite nonzero direction.
    Ray(PointView Origin, PointView Direction);

    /// @brief Returns whether this positive ray intersects a closed axis-aligned box.
    /// @param rLowerBound Lower box bound.
    /// @param rUpperBound Upper box bound.
    /// @return True when the positive ray intersects the box.
    [[nodiscard]] bool IntersectsAabb(PointView rLowerBound, PointView rUpperBound) const noexcept;

    /// @brief Intersects this positive ray with one triangle.
    /// @param Triangle Triangle geometry.
    /// @param MinimumRayDistance Minimum physical ray distance treated as non-negative.
    /// @return Named intersection status and hit data.
    [[nodiscard]] RayTriangleIntersection
        IntersectTriangle(TriangleProxy<WithoutNormals> Triangle, double MinimumRayDistance) const noexcept;

    /// @brief Returns whether the ray direction is parallel to the triangle plane within an angular tolerance.
    /// @param Triangle Triangle geometry.
    /// @param Tolerance Non-negative dimensionless angular tolerance.
    /// @return True for parallel or exactly degenerate triangles.
    [[nodiscard]] bool IsParallel(TriangleProxy<WithoutNormals> Triangle) const noexcept(NOTDEBUG);

private:
    /// @brief Intersects an AABB for directions containing zero or negative components.
    [[nodiscard]] bool IntersectsGeneralAabb(PointView rLowerBound, PointView rUpperBound) const noexcept;

    /// @brief Intersects an AABB using the strictly-positive-direction fast path.
    [[nodiscard]] bool IntersectsPositiveAabb(PointView rLowerBound, PointView rUpperBound) const noexcept;

    PointType mOrigin{};  ///< Finite ray origin.
    Vector3d mDirection{};  ///< Normalized ray direction.
    Vector3d mInverseDirection{};  ///< Reciprocal direction used by slab tests.
    Vector3i mDirectionSigns{};  ///< Lower/upper bound selector for each direction sign.
    bool mDirectionIsStrictlyPositive{};  ///< Whether every direction component is positive.
};

}  // namespace queso::embedding

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

//// External includes
#include "aabb_tree/AABB_base.h"

//// Project includes
#include "queso/containers/triangle_proxies.hpp"
#include "queso/includes/define.hpp"

namespace queso::embedding::detail {

/// @brief Exact separating-axis predicate between one axis-aligned box and one triangle.
/// @details Candidate generation remains conservative; this predicate decides exact intersection. With zero inset,
///          touching counts as intersection. A positive inset erodes every box half-extent and excludes touching-only
///          geometry.
class AabbTriangleIntersector : private aabb_base::AABB_base
{
public:
    /// @brief Constructs an intersector for one axis-aligned box.
    /// @param rLowerBound Lower box bound.
    /// @param rUpperBound Upper box bound.
    AabbTriangleIntersector(PointView rLowerBound, PointView rUpperBound)
        : aabb_base::AABB_base(
              Vector3d{ rLowerBound[0], rLowerBound[1], rLowerBound[2] },
              Vector3d{ rUpperBound[0], rUpperBound[1], rUpperBound[2] }
          )
    {}

    /// @brief Returns whether the inset box intersects a triangle.
    /// @param Triangle Triangle geometry.
    /// @param Inset Distance subtracted from every box half-extent.
    /// @return True when the inset box and triangle intersect.
    [[nodiscard]] bool IntersectsTriangle(TriangleProxy<WithoutNormals> Triangle, double Inset) const noexcept;

private:
    /// @brief Tests overlap after projecting the box and triangle onto one separating-axis candidate.
    /// @param rAxisX Unit x axis.
    /// @param rAxisY Unit y axis.
    /// @param rAxisZ Unit z axis.
    /// @param rFirst First translated triangle vertex.
    /// @param rSecond Second translated triangle vertex.
    /// @param rThird Third translated triangle vertex.
    /// @param rExtent Box half-extents after inset.
    /// @param rTestAxis Separating-axis candidate.
    /// @return True when projected intervals overlap.
    [[nodiscard]] bool CheckAxis(
        PointView rAxisX,
        PointView rAxisY,
        PointView rAxisZ,
        PointView rFirst,
        PointView rSecond,
        PointView rThird,
        PointView rExtent,
        PointView rTestAxis
    ) const noexcept;
};

}  // namespace queso::embedding::detail

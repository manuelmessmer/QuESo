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
#include <algorithm>
#include <optional>

//// Project includes
#include "queso/containers/triangle_mesh_view.hpp"
#include "queso/includes/define.hpp"
#include "queso/utilities/mapping_utilities.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class CellDomain
/// @author Manuel Messmer
/// @brief Active-domain representation for untrimmed elements.
/// TODO: Add domain concept.
class CellDomain
{
public:
    ///@name Life cycle
    ///@{

    /// Default constructor.
    CellDomain() = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns active-domain bounds in the requested coordinate space.
    /// @tparam TSpace Coordinate space.
    /// @param rBounds Element-local bounds.
    /// @return BoundingBoxType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] const BoundingBoxType& GetBounds(const ElementBounds& rBounds) const noexcept
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return rBounds.global;
        } else {
            return rBounds.parametric;
        }
    }

    /// @brief Returns the boundary mesh of the full-cell domain.
    /// @return TriangleMeshView
    [[nodiscard]] TriangleMeshView GetBoundaryMesh(const ElementBounds& rBounds) const
    { return GetOrCreateBoundaryMesh(rBounds).View(); }

    /// @brief Returns boundary integration points of the full-cell domain in the coordinate space specified by TSpace.
    /// @tparam TBoundaryIntegrationPointType
    /// @tparam TSpace Coordinate space.
    /// @return Vector of boundary integration points.
    template<typename TBoundaryIntegrationPointType, CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] std::vector<TBoundaryIntegrationPointType> GetBoundaryIps(const ElementBounds& rBounds) const
    {
        std::vector<TBoundaryIntegrationPointType> boundary_ips{};
        const auto boundary_mesh = GetOrCreateBoundaryMesh(rBounds).View();

        boundary_ips.reserve(boundary_mesh.NumOfTriangles() * 12UL);
        constexpr IndexType method = 3;
        boundary_mesh.template VisitEachTriangle<WithNormals>([&](const auto& rTriangle) {
            auto new_points = TriangleUtilities::GetIPsGlobal<TBoundaryIntegrationPointType>(rTriangle, method);
            if constexpr (TSpace == CoordinateSpace::global) {
                std::ranges::copy(new_points, std::back_inserter(boundary_ips));
            } else {
                std::ranges::transform(new_points, std::back_inserter(boundary_ips), [&](const auto& rPoint) {
                    return mapping::ToParametric(rPoint, rBounds);
                });
            }
        });

        return boundary_ips;
    }

    /// @brief Returns true if the query point is inside the domain bounds.
    /// @tparam TSpace Coordinate space.
    /// @param rPoint Query point in the coordinate space specified by TSpace.
    /// @return bool
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] bool IsInside(const PointType& rPoint, const ElementBounds& rBounds) const noexcept
    {
        const auto& bounds = GetBounds<TSpace>(rBounds);
        return bounds.lower[0] <= rPoint[0] && rPoint[0] <= bounds.upper[0] && bounds.lower[1] <= rPoint[1]
               && rPoint[1] <= bounds.upper[1] && bounds.lower[2] <= rPoint[2] && rPoint[2] <= bounds.upper[2];
    }

    /// @brief Returns intersection state of a query AABB against the full-cell domain.
    /// @tparam TSpace Coordinate space.
    /// @param rLowerBound Lower bound of query AABB in the coordinate space specified by TSpace.
    /// @param rUpperBound Upper bound of query AABB in the coordinate space specified by TSpace.
    /// @param Tolerance Query-box shrink tolerance.
    /// @return IntersectionStateType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        const ElementBounds& rBounds,
        double Tolerance = SNAPTOL
    ) const noexcept
    {
        const double lower_x = rLowerBound[0] + Tolerance;
        const double lower_y = rLowerBound[1] + Tolerance;
        const double lower_z = rLowerBound[2] + Tolerance;
        const double upper_x = rUpperBound[0] - Tolerance;
        const double upper_y = rUpperBound[1] - Tolerance;
        const double upper_z = rUpperBound[2] - Tolerance;

        const auto& bounds = GetBounds<TSpace>(rBounds);
        const bool disjoint = upper_x < bounds.lower[0] || bounds.upper[0] < lower_x || upper_y < bounds.lower[1]
                              || bounds.upper[1] < lower_y || upper_z < bounds.lower[2] || bounds.upper[2] < lower_z;
        if (disjoint) { return IntersectionState::outside; }

        const bool inside = bounds.lower[0] <= lower_x && upper_x <= bounds.upper[0] && bounds.lower[1] <= lower_y
                            && upper_y <= bounds.upper[1] && bounds.lower[2] <= lower_z && upper_z <= bounds.upper[2];
        return inside ? IntersectionState::inside : IntersectionState::trimmed;
    }

    ///@}
private:
    ///@name Private operations
    ///@{

    /// @brief Returns lazily initialized boundary mesh of the full cell.
    /// @return const TriangleMesh&
    [[nodiscard]] const TriangleMesh& GetOrCreateBoundaryMesh(const ElementBounds& rBounds) const
    {
        if (!mBoundaryMesh) {
            const auto& bounds = rBounds.global;
            mBoundaryMesh.emplace(MeshUtilities::MakeMeshBox(bounds.lower, bounds.upper));
        }
        return *mBoundaryMesh;
    }
    ///@}
    ///@name Private members
    ///@{

    mutable std::optional<TriangleMesh> mBoundaryMesh{};
    ///@}
};
///@}
}// namespace queso

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

/// STL includes
#include <algorithm>
#include <iterator>
/// Project includes
#include "queso/embedding/geometry_query.h"
#include "queso/embedding/trimmed_domain_on_plane.h"
#include "queso/utilities/mapping_utilities.hpp"
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// @class  TrimmedDomain
/// @author Manuel Messmer
/// @brief  Provides geometric operations for a trimmed domain, such as constructing a closed boundary mesh from a
///         clipped triangle mesh.
/// @details Uses an AABB tree for fast spatial queries.
/// TODO: Add domain concept.
class TrimmedDomain
{

public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor for a trimmed domain.
    /// @param rTriangleMesh Clipped triangle mesh. It must contain edges on the boundary planes of the enclosing AABB
    ///                      (see TrimmedDomainOnPlane).
    /// @param rLowerBound Lower bound of the enclosing AABB.
    /// @param rUpperBound Upper bound of the enclosing AABB.
    /// @param pOperator Pointer to the BRepOperator used as a fallback for IsInside() queries if
    ///                  IsInsideTrimmedDomain() is inconclusive.
    /// @param MinNumberOfTriangles Minimum number of triangles used to discretize the boundary of this trimmed domain.
    /// @param SwitchPlaneOrientation If true, the orientation of edges on TrimmedDomainOnPlane is switched.
    TrimmedDomain(
        ClippedTriangleMesh&& rTriangleMesh,
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        const BRepOperator* pOperator,
        IndexType MinNumberOfTriangles = 100,
        bool SwitchPlaneOrientation = false
    )
        : mpBrepOperatorGlobal(pOperator), mClippedMesh(std::move(rTriangleMesh)), mClosedMesh(mClippedMesh.Mesh()),
          mGeometryQuery(mClippedMesh.MeshView(), false)
    {
        // Set relative snap tolerance.
        mSnapTolerance = RelativeSnapTolerance(rLowerBound, rUpperBound);

        // Construct trimmed domain on plane upper bound of AABB.
        bool upper_bound = true;
        auto p_trimmed_domain_upper_x =
            MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_upper_y =
            MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_upper_z =
            MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);
        // Construct trimmed domain on plane lower bound of AABB.
        upper_bound = false;
        auto p_trimmed_domain_lower_x =
            MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_lower_y =
            MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_lower_z =
            MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, rLowerBound, rUpperBound, this, SwitchPlaneOrientation);

        if (mClippedMesh.NumOfTriangles() != 0) {
            auto p_t1 = p_trimmed_domain_lower_x->pGetTriangulation(mClippedMesh);
            auto p_t2 = p_trimmed_domain_upper_x->pGetTriangulation(mClippedMesh);
            auto p_t3 = p_trimmed_domain_lower_y->pGetTriangulation(mClippedMesh);
            auto p_t4 = p_trimmed_domain_upper_y->pGetTriangulation(mClippedMesh);
            auto p_t5 = p_trimmed_domain_lower_z->pGetTriangulation(mClippedMesh);
            auto p_t6 = p_trimmed_domain_upper_z->pGetTriangulation(mClippedMesh);

            const IndexType num_triangles = p_t1->NumOfTriangles() + p_t2->NumOfTriangles() + p_t3->NumOfTriangles()
                                            + p_t4->NumOfTriangles() + p_t5->NumOfTriangles() + p_t6->NumOfTriangles();

            mClosedMesh.Reserve(2UL * num_triangles);

            MeshUtilities::Append(mClosedMesh, p_t1->Mesh());
            MeshUtilities::Append(mClosedMesh, p_t2->Mesh());
            MeshUtilities::Append(mClosedMesh, p_t3->Mesh());
            MeshUtilities::Append(mClosedMesh, p_t4->Mesh());
            MeshUtilities::Append(mClosedMesh, p_t5->Mesh());
            MeshUtilities::Append(mClosedMesh, p_t6->Mesh());

            MeshUtilities::Refine(mClosedMesh, MinNumberOfTriangles);
        }
    }

    /// Destructor
    ~TrimmedDomain() = default;
    /// Copy constructor
    TrimmedDomain(const TrimmedDomain& rOther) = delete;
    /// Assignment operator
    TrimmedDomain& operator=(const TrimmedDomain& rOther) = delete;

    /// Move constructor
    TrimmedDomain(TrimmedDomain&& rOther) noexcept
        : mpBrepOperatorGlobal(rOther.mpBrepOperatorGlobal), mClippedMesh(std::move(rOther.mClippedMesh)),
          mClosedMesh(std::move(rOther.mClosedMesh)), mGeometryQuery(mClippedMesh.MeshView(), false),
          mSnapTolerance(rOther.mSnapTolerance)
    {}

    /// Move assignment operator
    TrimmedDomain& operator=(TrimmedDomain&& rOther) noexcept
    {
        if (this != &rOther) {
            mpBrepOperatorGlobal = rOther.mpBrepOperatorGlobal;
            mClippedMesh = std::move(rOther.mClippedMesh);
            mClosedMesh = std::move(rOther.mClosedMesh);
            mGeometryQuery = GeometryQuery(mClippedMesh.MeshView(), false);
            mSnapTolerance = rOther.mSnapTolerance;
        }
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns a view of the closed boundary mesh of the trimmed domain.
    /// @return TriangleMeshView
    [[nodiscard]] TriangleMeshView GetBoundaryMesh() const
    { return mClosedMesh.View(); }

    /// @brief Returns boundary integration points of the active domain.
    /// @tparam BoundaryIntegrationPointType Boundary integration-point type.
    /// @tparam TSpace Coordinate space of the returned integration points.
    /// @param rBounds Element-local bounds used for coordinate transformations.
    /// @return std::vector<BoundaryIntegrationPointType>
    template<typename BoundaryIntegrationPointType, CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] std::vector<BoundaryIntegrationPointType>
        GetBoundaryIps([[maybe_unused]] const ElementBounds& rBounds) const
    {
        std::vector<BoundaryIntegrationPointType> result{};
        result.reserve(mClosedMesh.NumOfTriangles() * 12UL);

        mClosedMesh.View().VisitEachTriangle<WithNormals>([&](const auto& triangle) {
            constexpr IndexType method = 3;

            auto points = TriangleUtilities::GetIPsGlobal<BoundaryIntegrationPointType>(triangle, method);

            if constexpr (TSpace == CoordinateSpace::global) {
                result.insert(
                    result.end(), std::make_move_iterator(points.begin()), std::make_move_iterator(points.end())
                );
            } else {
                std::ranges::transform(points, std::back_inserter(result), [&](auto& rPoint) {
                    return mapping::ToParametric(rPoint, rBounds);
                });
            }
        });

        return result;
    }

    /// @brief Returns the trimmed-domain bounds in the requested coordinate space.
    /// @tparam TSpace Coordinate space.
    /// @param rBounds Element-local bounds used for coordinate transformations.
    /// @return BoundingBoxType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] BoundingBoxType GetBounds([[maybe_unused]] const ElementBounds& rBounds) const
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return GetBoundingBoxOfTrimmedDomain();
        } else {
            const auto global_bounds = GetBoundingBoxOfTrimmedDomain();
            return MakeBox(
                mapping::ToParametric(global_bounds.lower, rBounds), mapping::ToParametric(global_bounds.upper, rBounds)
            );
        }
    }

    /// @brief Returns whether a global-space point lies inside the trimmed domain.
    /// @details Convenience overload for global-space queries that avoids passing
    /// element bounds.
    /// @tparam TSpace Must be CoordinateSpace::global.
    /// @param rPoint Query point in global coordinates.
    /// @return bool
    template<CoordinateSpace TSpace = CoordinateSpace::global>
        requires(TSpace == CoordinateSpace::global)
    [[nodiscard]] bool IsInside(const PointType& rPoint) const
    { return IsInsideTrimmedDomain(rPoint); }

    /// @brief Returns whether a point lies inside the trimmed domain.
    /// @tparam TSpace Coordinate space of the query point.
    /// @param rPoint Query point.
    /// @param rBounds Element-local bounds used for coordinate transformations.
    /// @return bool
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] bool IsInside(const PointType& rPoint, [[maybe_unused]] const ElementBounds& rBounds) const
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return IsInsideTrimmedDomain(rPoint);
        } else {
            return IsInsideTrimmedDomain(mapping::ToGlobal(rPoint, rBounds));
        }
    }

    /// @brief Returns the intersection state of a query box against the trimmed domain.
    /// @tparam TSpace Coordinate space of the query bounds.
    /// @param rLowerBound Lower bound of the query box.
    /// @param rUpperBound Upper bound of the query box.
    /// @param rBounds Element-local bounds used for coordinate transformations.
    /// @param Tolerance Shrink tolerance applied to the query box.
    /// @return IntersectionStateType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        const ElementBounds& rBounds,
        double Tolerance = SNAPTOL
    ) const
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return GetIntersectionState(rLowerBound, rUpperBound, Tolerance);
        } else {
            return GetIntersectionState(
                mapping::ToGlobal(rLowerBound, rBounds), mapping::ToGlobal(rUpperBound, rBounds), Tolerance
            );
        }
    }

    ///@}
private:
    ///@}
    ///@name Private Operations
    ///@{

    /// @brief Returns the bounding box of the trimmed domain.
    /// @details The result may be smaller than the full element domain.
    /// @return BoundingBoxType
    [[nodiscard]] BoundingBoxType GetBoundingBoxOfTrimmedDomain() const;

    /// @brief Returns whether a global-space point lies inside the trimmed domain.
    /// @details Performs a fast local test first via
    ///          IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess). If that test is inconclusive,
    ///          a fallback query is performed using the BRepOperator when available.
    /// @param rPoint Query point in global coordinates.
    /// @return bool
    [[nodiscard]] bool IsInsideTrimmedDomain(const PointType& rPoint) const;

    /// @brief Returns whether a global-space point lies inside the trimmed domain.
    /// @details Expects the point to lie inside the enclosing AABB; that check is omitted here. The test performs
    ///          ray tracing in the direction of the first triangle and searches for all intersections of the ray.
    ///          Inside/outside is determined from the orientation of the closest intersected triangle.
    /// @param rPoint Query point in global coordinates.
    /// @param[out] rSuccess Set to false if the result is ambiguous, for example when all candidate triangles are
    ///                      parallel to the tracing ray.
    /// @return bool
    [[nodiscard]] bool IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess) const;

    /// @brief Returns the intersection state of an AABB against the trimmed domain.
    /// @note This test is performed only on `mClippedMesh` for efficiency. This is primarily used by the octree.
    /// @param rLowerBound Lower bound of the query AABB.
    /// @param rUpperBound Upper bound of the query AABB.
    /// @param Tolerance Shrink tolerance applied to the query AABB. If `Tolerance == 0`, touching counts as an
    ///                  intersection. If `Tolerance > 0`, touching is not treated as an intersection.
    /// @return IntersectionStateType with values `Inside`, `Outside`, or `Trimmed`.
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        double Tolerance = SNAPTOL
    ) const;

    ///@}
    ///@name Private Members
    ///@{

    const BRepOperator* mpBrepOperatorGlobal;

    ClippedTriangleMesh mClippedMesh;
    TriangleMesh mClosedMesh;
    GeometryQuery mGeometryQuery;
    double mSnapTolerance;
    ///@}
};
///@}
}// namespace queso

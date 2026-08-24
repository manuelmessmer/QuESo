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
#include <cmath>
#include <functional>
#include <iterator>
#include <span>
#include <vector>
/// Project includes
#include "queso/embedding/cell_face_closure.h"
#include "queso/embedding/local_surface_classifier.h"
#include "queso/embedding/mesh_query.h"
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

    /// @brief Constructs a trimmed domain with a construction-only global fallback classifier.
    /// @details The fallback is called only when a local clipped-surface query is inconclusive. It is not retained by
    ///          the completed domain.
    /// @param rSection Source surface, exact cell bounds, contours, and tolerance.
    /// @param rGlobalIsInside Closed-source-mesh point classifier available during construction.
    /// @param MinNumberOfTriangles Minimum number of triangles used to discretize the completed boundary.
    TrimmedDomain(
        embedding::CellSurfaceSection&& rSection,
        const std::function<bool(PointView)>& rGlobalIsInside,
        IndexType MinNumberOfTriangles = 100
    )
        : TrimmedDomain(PrepareClosure(std::move(rSection), rGlobalIsInside), MinNumberOfTriangles)
    {}

    /// Destructor
    ~TrimmedDomain() = default;
    /// Copy constructor
    TrimmedDomain(const TrimmedDomain& rOther) = delete;
    /// Assignment operator
    TrimmedDomain& operator=(const TrimmedDomain& rOther) = delete;

    /// Move constructor
    TrimmedDomain(TrimmedDomain&& rOther)
        : mClippedMesh(std::move(rOther.mClippedMesh)), mClosedMesh(std::move(rOther.mClosedMesh)),
          mMeshQuery(mClippedMesh.View(), embedding::MeshQueryMode::OrientedSurface, rOther.mGeometryTolerance),
          mGeometryTolerance(rOther.mGeometryTolerance), mActiveBounds(rOther.mActiveBounds)
    {}

    /// Move assignment operator
    TrimmedDomain& operator=(TrimmedDomain&& rOther)
    {
        if (this != &rOther) {
            mClippedMesh = std::move(rOther.mClippedMesh);
            mClosedMesh = std::move(rOther.mClosedMesh);
            mMeshQuery = embedding::MeshQuery(
                mClippedMesh.View(), embedding::MeshQueryMode::OrientedSurface, rOther.mGeometryTolerance
            );
            mGeometryTolerance = rOther.mGeometryTolerance;
            mActiveBounds = rOther.mActiveBounds;
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
            return mActiveBounds;
        } else {
            return MakeBox(
                mapping::ToParametric(mActiveBounds.lower, rBounds), mapping::ToParametric(mActiveBounds.upper, rBounds)
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
    /// @param Policy AABB boundary-touching policy.
    /// @return IntersectionStateType
    template<CoordinateSpace TSpace = CoordinateSpace::global>
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        const ElementBounds& rBounds,
        AabbIntersectionPolicy Policy = AabbIntersectionPolicy::Exact
    ) const
    {
        if constexpr (TSpace == CoordinateSpace::global) {
            return GetIntersectionState(rLowerBound, rUpperBound, Policy);
        } else {
            return GetIntersectionState(
                mapping::ToGlobal(rLowerBound, rBounds), mapping::ToGlobal(rUpperBound, rBounds), Policy
            );
        }
    }

    ///@}
private:
    ///@}
    ///@name Private Operations
    ///@{

    struct PreparedClosure
    {
        TriangleMesh clipped_surface;
        TriangleMesh closed_mesh;
        GeometryTolerance tolerance;
    };

    [[nodiscard]] static PreparedClosure
        PrepareClosure(embedding::CellSurfaceSection&& rSection, const std::function<bool(PointView)>& rGlobalIsInside)
    {
        const GeometryTolerance tolerance = rSection.GetGeometryTolerance();
        const BoundingBoxType bounds = rSection.CellBounds();
        TriangleMesh clipped_surface = rSection.TakeSurface();
        const auto IsInside = [&](PointView rPoint) {
            const PointType point{ rPoint[0], rPoint[1], rPoint[2] };
            const auto local =
                embedding::detail::ClassifyOnBoundedSide(point, clipped_surface.View(), tolerance.ZeroLength());
            if (local == embedding::detail::LocalSurfaceClassification::Inside) { return true; }
            if (local == embedding::detail::LocalSurfaceClassification::Outside) { return false; }
            // The global classifier exists only during construction and is never retained by the completed domain.
            return rGlobalIsInside ? rGlobalIsInside(point) : false;
        };
        const auto BuildClosedMesh = [&](bool SwitchAxes) {
            TriangleMesh result = clipped_surface;
            for (IndexType face_index = 0; face_index < 6; ++face_index) {
                const embedding::CellFace face = embedding::CellFace::FromIndex(face_index);
                const TriangleMesh closure = embedding::detail::BuildCellFaceClosure(
                    rSection.FaceSegments(face), bounds, face, tolerance, IsInside, SwitchAxes
                );
                MeshUtilities::Append(result, closure);
            }
            return result;
        };

        TriangleMesh closed_mesh = BuildClosedMesh(false);
        const double fixed_frame_quality = MeshUtilities::EstimateQuality(closed_mesh.View());
        TriangleMesh switched_mesh = BuildClosedMesh(true);
        const double swapped_frame_quality = MeshUtilities::EstimateQuality(switched_mesh.View());
        if (swapped_frame_quality < fixed_frame_quality) { closed_mesh = std::move(switched_mesh); }
        return { std::move(clipped_surface), std::move(closed_mesh), tolerance };
    }

    TrimmedDomain(PreparedClosure&& rPrepared, IndexType MinNumberOfTriangles)
        : mClippedMesh(std::move(rPrepared.clipped_surface)), mClosedMesh(std::move(rPrepared.closed_mesh)),
          mMeshQuery(mClippedMesh.View(), embedding::MeshQueryMode::OrientedSurface, rPrepared.tolerance),
          mGeometryTolerance(rPrepared.tolerance)
    {
        MeshUtilities::Refine(mClosedMesh, MinNumberOfTriangles, mGeometryTolerance.ZeroArea());
        mActiveBounds = ComputeBounds(mClosedMesh);
    }

    [[nodiscard]] static BoundingBoxType ComputeBounds(const TriangleMesh& rMesh)
    {
        const auto [lower, upper] = MeshUtilities::BoundingBox(rMesh.View());
        return MakeBox(lower, upper);
    }

    /// @brief Returns whether a global-space point lies inside the trimmed domain.
    /// @details Performs a deterministic local clipped-section query. Inconclusive queries classify as outside.
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
    /// @param Policy AABB boundary-touching policy.
    /// @return IntersectionStateType with values `Inside`, `Outside`, or `Trimmed`.
    [[nodiscard]] IntersectionStateType GetIntersectionState(
        const PointType& rLowerBound,
        const PointType& rUpperBound,
        AabbIntersectionPolicy Policy = AabbIntersectionPolicy::Exact
    ) const;

    ///@}
    ///@name Private Members
    ///@{

    TriangleMesh mClippedMesh;
    TriangleMesh mClosedMesh;
    embedding::MeshQuery mMeshQuery;
    GeometryTolerance mGeometryTolerance;
    BoundingBoxType mActiveBounds;
    ///@}
};
///@}
}  // namespace queso

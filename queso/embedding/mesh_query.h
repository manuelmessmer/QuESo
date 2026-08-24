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
#include <optional>
#include <span>
#include <vector>

//// Project includes
#include "queso/containers/geometry_tolerance.hpp"
#include "queso/containers/triangle_mesh_view.hpp"
#include "queso/embedding/aabb_tree.h"
#include "queso/embedding/ray.h"
#include "queso/includes/aabb_intersection_policy.hpp"

namespace queso::embedding {

/// @brief Ray-classification policy fixed for one MeshQuery.
enum class MeshQueryMode : std::uint8_t {
    Closed,  ///< Classify by odd/even intersection parity of a closed surface.
    OrientedSurface  ///< Classify from the winding of the nearest intersected open surface.
};

/// @brief Result of classifying one ray origin against a mesh.
struct RayClassification
{
    bool is_inside;  ///< Whether the conclusive classification is inside or on the bounded side.
    bool is_conclusive;  ///< Whether the ray avoided ambiguous edge, vertex, and parallel configurations.
};

/// @class MeshQuery
/// @brief Provides immutable spatial queries against one triangle mesh.
/// @details Owns the non-owning mesh view and a snapshot AABB acceleration structure. The referenced mesh must outlive
///          the query, and its geometry, topology, triangle order, and triangle count must remain unchanged. Concurrent
///          const queries require thread-safe const access to the underlying mesh. Candidate IDs are valid only for the
///          query that produced them.
class MeshQuery
{
public:
    ///@name Life Cycle
    ///@{

    /// @brief Constructs the query and its acceleration structure.
    /// @param rMesh Triangle mesh view. The referenced mesh must outlive this query.
    /// @param Mode Fixed ray-classification policy.
    /// @param rTolerance Immutable geometry tolerance for this query's coordinate scale.
    MeshQuery(const TriangleMeshView& rMesh, MeshQueryMode Mode, const GeometryTolerance& rTolerance)
        : mTriangleMesh(rMesh), mTree(rMesh), mMode(Mode), mGeometryTolerance(rTolerance)
    {}

    /// @brief MeshQuery is non-copyable because its acceleration structure is an expensive mesh-specific snapshot.
    MeshQuery(const MeshQuery&) = delete;
    /// @brief MeshQuery is non-copyable because its mesh view and acceleration structure must remain paired.
    MeshQuery& operator=(const MeshQuery&) = delete;
    /// @brief Transfers a query and its acceleration structure.
    MeshQuery(MeshQuery&&) = default;
    /// @brief Transfers a query and its acceleration structure.
    /// @return Reference to this query.
    MeshQuery& operator=(MeshQuery&&) = default;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Returns whether a point lies within the global mesh bounding box.
    /// @param rPoint Query point.
    /// @return True when the point lies within the bounding box.
    [[nodiscard]] bool IsWithinBoundingBox(PointView rPoint) const noexcept;

    /// @brief Returns the global mesh bounding box when the mesh is non-empty.
    /// @return Bounding box of all triangles, or std::nullopt for an empty mesh.
    [[nodiscard]] std::optional<BoundingBoxType> BoundingBox() const noexcept;

    /// @brief Returns conservative triangle candidates for an axis-aligned bounding box.
    /// @details Bounds must be finite and non-inverted. Debug builds validate this precondition. Candidate order is
    ///          unspecified.
    /// @param rLowerBound Lower AABB bound.
    /// @param rUpperBound Upper AABB bound.
    /// @return Candidate triangle IDs. Candidates require an exact intersection test before use as intersections.
    [[nodiscard]] std::vector<IndexType> GetAabbCandidates(PointView rLowerBound, PointView rUpperBound) const;

    /// @brief Exact-tests supplied triangle candidates against an axis-aligned bounding box.
    /// @details This operation does not query the AABB tree. Bounds must be finite and non-inverted; snap-eroded
    ///          queries must retain positive extent. Debug builds validate these preconditions. Exact is the default
    ///          and counts touching geometry as an intersection. Snap-eroded policy excludes geometry that only
    ///          touches the query-box boundary.
    /// @param rCandidateIds Candidate triangle IDs produced by this query.
    /// @param rLowerBound Lower AABB bound.
    /// @param rUpperBound Upper AABB bound.
    /// @param Policy AABB boundary-touching policy.
    /// @return True when any supplied candidate intersects the AABB.
    [[nodiscard]] bool IntersectsAabb(
        std::span<const IndexType> rCandidateIds,
        PointView rLowerBound,
        PointView rUpperBound,
        AabbIntersectionPolicy Policy = AabbIntersectionPolicy::Exact
    ) const;

    /// @brief Returns supplied candidate IDs that exactly intersect an axis-aligned bounding box.
    /// @details This operation does not query the AABB tree. Bounds must be finite and non-inverted; debug builds
    ///          validate this precondition. Touching geometry counts as an intersection.
    /// @param rCandidateIds Candidate triangle IDs produced by this query.
    /// @param rLowerBound Lower AABB bound.
    /// @param rUpperBound Upper AABB bound.
    /// @return Exactly intersecting triangle IDs.
    [[nodiscard]] std::vector<IndexType> GetIntersectedTriangleIds(
        std::span<const IndexType> rCandidateIds,
        PointView rLowerBound,
        PointView rUpperBound
    ) const;

    /// @brief Classifies a ray origin according to the mode selected at construction.
    /// @details Closed mode uses positive-ray parity and ignores exactly degenerate triangles. Oriented-surface mode
    ///          uses the nearest hit's winding. A near-origin hit is conclusively outside; edge, vertex, and parallel
    ///          configurations are inconclusive; a no-hit ray is conclusively outside.
    /// @param rRay Query ray.
    /// @return Named classification result, which also supports structured binding.
    [[nodiscard]] RayClassification Classify(const Ray& rRay) const;

    /// @brief Returns the non-owning mesh view used by this query.
    /// @return Triangle mesh view.
    [[nodiscard]] TriangleMeshView MeshView() const noexcept
    { return mTriangleMesh; }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /// @brief Classifies by odd/even positive-ray intersection parity.
    /// @param rRay Valid normalized ray.
    /// @return Closed-surface classification and conclusiveness.
    [[nodiscard]] RayClassification ClassifyClosed(const Ray& rRay) const;

    /// @brief Classifies from the winding of the nearest positive-ray intersection.
    /// @param rRay Valid normalized ray.
    /// @return Oriented-surface classification and conclusiveness.
    [[nodiscard]] RayClassification ClassifyOrientedSurface(const Ray& rRay) const;

    ///@}
    ///@name Private Members
    ///@{

    TriangleMeshView mTriangleMesh;  ///< Non-owning immutable source mesh view.
    detail::AabbTree mTree;  ///< Mesh-specific acceleration snapshot.
    MeshQueryMode mMode;  ///< Fixed ray-classification policy.
    GeometryTolerance mGeometryTolerance;  ///< Immutable physical tolerance for this query's coordinate scale.

    ///@}
};

}  // namespace queso::embedding

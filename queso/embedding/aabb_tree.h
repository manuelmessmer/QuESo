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
#include <vector>

//// External includes
#include "aabb_tree/AABB_base.h"

//// Project includes
#include "queso/containers/triangle_mesh_view.hpp"
#include "queso/includes/define.hpp"

namespace queso::embedding {

class Ray;

namespace detail {

    /// @brief Immutable AABB acceleration structure over one triangle mesh view.
    /// @details Adapts the external dynamic tree for construction-once, query-only use. Candidate queries are
    ///          conservative and return triangle IDs in tree-traversal order. The referenced mesh must remain
    ///          unchanged and outlive the tree.
    class AabbTree : private aabb_base::Tree_base
    {
    public:
        /// @brief Builds one tree node per source triangle bounding box.
        /// @param rTriangleMesh Immutable source mesh view.
        explicit AabbTree(const TriangleMeshView& rTriangleMesh);

        /// @brief Returns whether a point lies inside the complete mesh bounding box.
        /// @param rPoint Query point.
        /// @return True when the point lies within the closed bounds.
        [[nodiscard]] bool IsWithinBoundingBox(PointView rPoint) const noexcept;

        /// @brief Returns the complete triangle-mesh bounding box.
        /// @details An empty tree retains sentinel bounds; MeshQuery exposes empty bounds as std::nullopt.
        /// @return Mesh bounds or empty-tree sentinels.
        [[nodiscard]] BoundingBoxType BoundingBox() const noexcept
        { return { mLowerBound, mUpperBound }; }

        /// @brief Returns conservative candidates overlapping an axis-aligned query box.
        /// @param rLowerBound Lower query bound.
        /// @param rUpperBound Upper query bound.
        /// @return Candidate triangle IDs in tree-traversal order.
        [[nodiscard]] std::vector<IndexType> Query(PointView rLowerBound, PointView rUpperBound) const;

        /// @brief Returns conservative candidates intersected by a ray.
        /// @param rRay Ray query.
        /// @return Candidate triangle IDs in tree-traversal order.
        [[nodiscard]] std::vector<IndexType> Query(const Ray& rRay) const;

    private:
        using BaseTreeType = aabb_base::Tree_base;

        /// @brief Traverses tree nodes accepted by rIntersects and collects leaf triangle IDs.
        /// @tparam TPredicate Callable accepting lower and upper node bounds.
        /// @param rIntersects Node-intersection predicate.
        /// @return Candidate triangle IDs in deterministic traversal order.
        template<typename TPredicate>
        [[nodiscard]] std::vector<IndexType> QueryImpl(TPredicate&& rIntersects) const;

        PointType mLowerBound{};  ///< Lower bound across all source triangles.
        PointType mUpperBound{};  ///< Upper bound across all source triangles.
    };

}  // namespace detail
}  // namespace queso::embedding

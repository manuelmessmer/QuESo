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

//// Own include
#include "queso/embedding/aabb_tree.h"

//// STL includes
#include <algorithm>

//// Project includes
#include "queso/embedding/ray.h"

namespace queso::embedding::detail {

AabbTree::AabbTree(const TriangleMeshView& rTriangleMesh) : aabb_base::Tree_base(3, 0.0, 16, false)
{
    mLowerBound = { MAXD, MAXD, MAXD };
    mUpperBound = { LOWESTD, LOWESTD, LOWESTD };
    rTriangleMesh.VisitEachTriangle<WithoutNormals>([&, TriangleId = 0U](const auto& rTriangle) mutable {
        const PointType x_values{ rTriangle.P1[0], rTriangle.P2[0], rTriangle.P3[0] };
        const PointType y_values{ rTriangle.P1[1], rTriangle.P2[1], rTriangle.P3[1] };
        const PointType z_values{ rTriangle.P1[2], rTriangle.P2[2], rTriangle.P3[2] };
        const auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
        const auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
        const auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());
        const PointType lower{ *x_min_max.first, *y_min_max.first, *z_min_max.first };
        const PointType upper{ *x_min_max.second, *y_min_max.second, *z_min_max.second };
        insertParticle(TriangleId++, lower, upper);
        for (IndexType axis = 0; axis < 3; ++axis) {
            mLowerBound[axis] = std::min(mLowerBound[axis], lower[axis]);
            mUpperBound[axis] = std::max(mUpperBound[axis], upper[axis]);
        }
    });
}

bool AabbTree::IsWithinBoundingBox(PointView rPoint) const noexcept
{
    if (rPoint[0] < mLowerBound[0] || rPoint[0] > mUpperBound[0] || rPoint[1] < mLowerBound[1]
        || rPoint[1] > mUpperBound[1] || rPoint[2] < mLowerBound[2] || rPoint[2] > mUpperBound[2]) {
        return false;
    }

    return true;
}


std::vector<IndexType> AabbTree::Query(PointView rLowerBound, PointView rUpperBound) const
{
    return QueryImpl([rLowerBound, rUpperBound](PointView rNodeLower, PointView rNodeUpper) {
        for (IndexType axis = 0; axis < 3; ++axis) {
            if (rNodeUpper[axis] < rLowerBound[axis] || rNodeLower[axis] > rUpperBound[axis]) { return false; }
        }
        return true;
    });
}

std::vector<IndexType> AabbTree::Query(const Ray& rRay) const
{
    return QueryImpl([&rRay](PointView rNodeLower, PointView rNodeUpper) {
        return rRay.IntersectsAabb(rNodeLower, rNodeUpper);
    });
}

template<typename TPredicate>
std::vector<IndexType> AabbTree::QueryImpl(TPredicate&& rIntersects) const
{
    std::vector<IndexType> stack;
    stack.reserve(256);
    stack.push_back(BaseTreeType::Root());

    std::vector<IndexType> particles;

    while (stack.size() > 0) {
        IndexType node = stack.back();
        stack.pop_back();

        if (node == NULL_NODE) continue;

        const auto& r_aabb_base = BaseTreeType::Nodes()[node].aabb_base;

        // Test for overlap between the AABBs.
        if (rIntersects(r_aabb_base.lowerBound, r_aabb_base.upperBound)) {
            // Check that we're at a leaf node.
            if (BaseTreeType::Nodes()[node].isLeaf()) {
                particles.push_back(BaseTreeType::Nodes()[node].particle);
            } else {
                stack.push_back(BaseTreeType::Nodes()[node].left);
                stack.push_back(BaseTreeType::Nodes()[node].right);
            }
        }
    }
    return particles;
}

}  // namespace queso::embedding::detail

// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// Project includes
#include "embedding/aabb_tree.h"

bool AABB_tree::IsWithinBoundingBox(const PointType& rPoint) const {
    if(   rPoint[0] < mLowerBound[0]
        || rPoint[0] > mUpperBound[0]
        || rPoint[1] < mLowerBound[1]
        || rPoint[1] > mUpperBound[1]
        || rPoint[2] < mLowerBound[2]
        || rPoint[2] > mUpperBound[2])
    {
        return false;
    }

    return true;
}


std::vector<IndexType> AABB_tree::Query(const AABB_primitive_base& rAABB_primitive) const
{
    std::vector<IndexType> stack;
    stack.reserve(256);
    stack.push_back(BaseTreeType::Root());

    std::vector<IndexType> particles;

    while (stack.size() > 0)
    {
        IndexType node = stack.back();
        stack.pop_back();

        // Copy the AABB_base and cast it to AABB_primitive.
        const auto& aabb = static_cast<const AABB_primitive&>(BaseTreeType::Nodes()[node].aabb_base);

        if (node == NULL_NODE) continue;

        // Test for overlap between the AABBs.
        if (rAABB_primitive.intersect(aabb) )
        {
            // Check that we're at a leaf node.
            if (BaseTreeType::Nodes()[node].isLeaf())
            {
                particles.push_back(BaseTreeType::Nodes()[node].particle);
            }
            else
            {
                stack.push_back(BaseTreeType::Nodes()[node].left);
                stack.push_back(BaseTreeType::Nodes()[node].right);
            }
        }
    }
    return particles;
}

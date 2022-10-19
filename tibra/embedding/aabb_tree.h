// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_tree_INCLUDE_H
#define AABB_tree_INCLUDE_H

/// External libraries
#include "aabb_tree/AABB_base.h"

/// Project includes
#include "embedding/aabb_primitive.h"
#include "geometries/triangle_mesh.h"


///@name TIBRA Classes
///@{

/**
 * @class  AABB_tree
 * @author Manuel Messmer
 * @brief  Derives from AABB_base, see: https://github.com/lohedges/aabbcc.
*/
class AABB_tree : public aabb_base::Tree_base {

public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef TriangleMesh::Vector3d PointType;
    typedef aabb_base::Tree_base BaseTreeType;
    typedef aabb_base::AABB_base BaseAABBType;
    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor.
    /// @param Dimension
    /// @param TriangleMesh
    /// @param SkinThickness The skin thickness for fattened AABBs, as a fraction
    ///                      of the AABB_base base length.
    /// @param TouchIsOverlap Does touching count as overlapping in query operations?
    AABB_tree(const TriangleMesh& TriangleMesh ) :
            aabb_base::Tree_base(3, 0.0, 16,  false)
    {
        const double max_limit = std::numeric_limits<double>::max();
        mLowerBound = {max_limit, max_limit, max_limit};
        mUpperBound = {-max_limit, -max_limit, -max_limit};
        for( int i = 0; i < TriangleMesh.NumOfTriangles(); ++i){
            const auto& p1 = TriangleMesh.P1(i);
            const auto& p2 = TriangleMesh.P2(i);
            const auto& p3 = TriangleMesh.P3(i);

            const PointType x_values = {p1[0], p2[0], p3[0]};
            const PointType y_values = {p1[1], p2[1], p3[1]};
            const PointType z_values = {p1[2], p2[2], p3[2]};

            auto [x_min, x_max] = std::minmax_element(x_values.begin(), x_values.end());
            auto [y_min, y_max] = std::minmax_element(y_values.begin(), y_values.end());
            auto [z_min, z_max] = std::minmax_element(z_values.begin(), z_values.end());

            std::array<double,3> lower_bound = {*x_min, *y_min, *z_min};
            std::array<double,3> upper_bound = {*x_max, *y_max, *z_max};
            this->insertParticle(i, lower_bound, upper_bound);

            mLowerBound[0] = std::min<double>(*x_min, mLowerBound[0]);
            mUpperBound[0] = std::max<double>(*x_max, mUpperBound[0]);

            mLowerBound[1] = std::min<double>(*y_min, mLowerBound[1]);
            mUpperBound[1] = std::max<double>(*y_max, mUpperBound[1]);

            mLowerBound[2] = std::min<double>(*z_min, mLowerBound[2]);
            mUpperBound[2] = std::max<double>(*z_max, mUpperBound[2]);
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Return true if point lies inside outer bounding box.
    ///@param rPoint
    ///@return bool
    bool IsWithinBoundingBox(const std::array<double, 3>& rPoint) const;

    ///@brief Get all potential interesections of Ray. Checks against bounding box of triangles.
    ///@param rRay
    ///@return std::vector<unsigned int> Holds Id's of triangles.
    ///@todo Should return ptr to std::vector<unsigned int>;
    std::vector<IndexType> Query(const AABB_primitive_base& rAABB_primitive) const;

    ///@}

private:
    ///@name Private Member variables
    ///@{
    std::array<double,3> mLowerBound{};
    std::array<double,3> mUpperBound{};
    ///@}
};

///@}

#endif //AABB_tree_INCLUDE_H
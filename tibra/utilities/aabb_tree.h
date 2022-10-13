// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_tree_INCLUDE_H
#define AABB_tree_INCLUDE_H

/// External includes
#include <memory>

/// External libraries
#include "aabb_tree/AABB_base.h"

/// Project includes
#include "geometries/triangle_mesh.h"


// Create alias for AABB
typedef aabb_base::AABB_base AABB;

constexpr double kEpsilon = 1e-14;

///@name TIBRA Classes
///@{


/**
 * @class  Ray
 * @author Manuel Messmer
 * @brief  Ray to be used in AABB_tree. Direction of ray must be positive oriented.
 *         Check is omitted for better performance.
*/
class Ray {
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef TriangleMesh::Vector3d Vector3d;
    typedef TriangleMesh::Vector3i Vector3i;

    ///@}
    ///@name Life cycle
    ///@{

    /// Constructor.
    Ray(const Vector3d& Origin, const Vector3d& Direction) :
        mOrigin(Origin), mDirection(Direction)
    {
        mInvDirection[0] = 1.0 / Direction[0];
        mInvDirection[1] = 1.0 / Direction[1];
        mInvDirection[2] = 1.0 / Direction[2];
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if origin of ray is inside aabb.
    ///@param aabb
    ///@return bool
    bool inside(AABB &aabb)
    {
        if(    mOrigin[0] < aabb.lowerBound[0]
            || mOrigin[1] < aabb.lowerBound[1]
            || mOrigin[2] < aabb.lowerBound[2]
            || mOrigin[0] > aabb.upperBound[0]
            || mOrigin[1] > aabb.upperBound[1]
            || mOrigin[2] > aabb.upperBound[2] ) {
                return false;
            }

        return true;
    }

    ///@brief Returns true if ray intersects aabb.
    ///@param aabb
    ///@param t distance to intersection.
    ///@return bool
    bool intersect(AABB &aabb, double &t)
    {
        double tmin, tmax, tymin, tymax, tzmin, tzmax;

        double lower_0 = aabb.lowerBound[0];
        double lower_1 = aabb.lowerBound[1];

        double upper_0 = aabb.upperBound[0];
        double upper_1 = aabb.upperBound[1];

        double origin_0 = mOrigin[0];
        double origin_1 = mOrigin[1];

        double inv_direction_0 = mInvDirection[0];
        double inv_direction_1 = mInvDirection[1];


        tmin = (lower_0 - origin_0) * inv_direction_0;
        tymax = (upper_1 - origin_1) * inv_direction_1;
        if(tmin > tymax) {
            return false;
        }

        tmax = (upper_0 - origin_0) * inv_direction_0;
        tymin = (lower_1 - origin_1) * inv_direction_1;


        if( tymin > tmax )
            return false;

        if (tymin > tmin)
            tmin = tymin;

        if (tymax < tmax)
            tmax = tymax;

        double lower_2 = aabb.lowerBound[2];
        double upper_2 = aabb.upperBound[2];
        double origin_2 = mOrigin[2];
        double inv_direction_2 = mInvDirection[2];


        tzmin = (lower_2 - origin_2) * inv_direction_2;
        if((tzmin > tmax))
            return false;

        tzmax = (upper_2 - origin_2) * inv_direction_2;
        if( tmin > tzmax )
            return false;

        if (tzmin > tmin)
            tmin = tzmin;

        t = tmin;
        if (t < 0) {
            if (tzmax < tmax)
                tmax = tzmax;
            t = tmax;
            return false;
        }

        return true;
    }

    ///@brief Returns true if ray intersects triangle (only checks intersections in positive direction).
    ///@param v0 Triangle Vertex 1
    ///@param v1 Triangle Vertex 2
    ///@param v2 Triangle Vertex 3
    ///@param t Distance to intersection.
    ///@param u Parametric coordinate 1.
    ///@param v Parametric coordinate 2.
    ///@return bool
    bool intersect(
        const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
        double &t, double &u, double &v)
    {
        // Substraction: v1-v0 and v2-v0
        Vector3d v0v1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        Vector3d v0v2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

        // Cross product: mDirection x v0v2
        Vector3d pvec = {mDirection[1]*v0v2[2] - mDirection[2]*v0v2[1],
                        mDirection[2]*v0v2[0] - mDirection[0]*v0v2[2],
                        mDirection[0]*v0v2[1] - mDirection[1]*v0v2[0]};

        // Dot product: v0v1 * pvec
        double det = v0v1[0]*pvec[0] + v0v1[1]*pvec[1] + v0v1[2]*pvec[2];

        if (fabs(det) < kEpsilon)
            return false;

        double invDet = 1 / det;

        // Substraction: mOrigin - v0
        Vector3d tvec = {mOrigin[0] - v0[0], mOrigin[1] - v0[1], mOrigin[2] - v0[2]};

        // Dot product x invDet: (tvec * pvec) * invDet
        u = (tvec[0]*pvec[0] + tvec[1]*pvec[1] + tvec[2]*pvec[2]) * invDet;

        if (u < 0 || u > 1)
            return false;

        // Cross product: tvec x v0v1
        Vector3d qvec = {tvec[1]*v0v1[2] - tvec[2]*v0v1[1],
                        tvec[2]*v0v1[0] - tvec[0]*v0v1[2],
                        tvec[0]*v0v1[1] - tvec[1]*v0v1[0]};

        // Dot product x invDet: (mDirection * qvec) * invDet
        v = (mDirection[0]*qvec[0] + mDirection[1]*qvec[1] + mDirection[2]*qvec[2]) * invDet;

        if (v < 0 || u + v > 1)
            return false;

        // Dot product: v0v2 * qvec
        t = (v0v2[0]*qvec[0] + v0v2[1]*qvec[1] + v0v2[2]*qvec[2]) * invDet;

        // Return true if ray intersects in negative direction.
        if( t < 0 )
            return false;

        return true;
    }

    ///@}
private:

    ///@name Private Members
    ///@{
    Vector3d mOrigin;
    Vector3d mDirection;
    Vector3d mInvDirection;
    Vector3i mSign;
    ///@}
};
///@} // End class Ray

///
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
            aabb_base::Tree_base(3, 0.000001, 16,  false)
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
    bool IsWithinBoundingBox(const std::array<double, 3>& rPoint){
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

    ///@brief Get all interesections of Ray
    ///@param rRay
    ///@return std::vector<unsigned int> Holds Id's of triangles.
    std::vector<unsigned int> Query(Ray& rRay)
    {
        std::vector<unsigned int> stack;
        stack.reserve(256);
        stack.push_back(BaseTreeType::Root());

        std::vector<unsigned int> particles;

        while (stack.size() > 0)
        {
            unsigned int node = stack.back();
            stack.pop_back();

            // Copy the AABB_base.
            AABB nodeAABB = BaseTreeType::Nodes()[node].aabb_base;

            if (node == NULL_NODE) continue;

            // Test for overlap between the AABBs.
            double t;
            if( rRay.inside(nodeAABB) ){
                // Check that we're at a leaf node.
                if (BaseTreeType::Nodes()[node].isLeaf() )
                {
                    particles.push_back(BaseTreeType::Nodes()[node].particle);
                }
                else
                {
                    stack.push_back(BaseTreeType::Nodes()[node].left);
                    stack.push_back(BaseTreeType::Nodes()[node].right);
                }
            }
            else if (rRay.intersect(nodeAABB, t) )
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
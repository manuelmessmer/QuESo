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
constexpr double kEpsilon = 1e-8;
///@name TIBRA Classes
///@{


/**
 * @class  Ray
 * @author Manuel Messmer
 * @brief  Ray to be used in AABB_tree
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
        mInvDirection[0] = 1 / Direction[0];
        mInvDirection[1] = 1 / Direction[1];
        mInvDirection[2] = 1 / Direction[2];

        mSign[0] = (mInvDirection[0] < 0);
        mSign[1] = (mInvDirection[1] < 0);
        mSign[2] = (mInvDirection[2] < 0);
    }

    ///@}
    ///@name Operations
    ///@{
    /// The code from this lesson returns intersections with the box that are either in front or
    /// behind the origin of the ray. For instance, if the ray's origin is inside the box (adjacent image),
    /// there will be two intersections: one in front of the ray and one behind. We know that an intersection
    /// is "behind" the origin of the ray when t is negative. When t is positive, the intersection is in front
    /// of the origin of the ray (with respect to the ray's direction).

    //bool intersect(const Ray &r, double &t) const
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

    bool intersect(AABB &aabb, double &t)
    {
        double tmin, tmax, tymin, tymax, tzmin, tzmax;
        //std::cout << "test " << std::endl;
        std::array<Vector3d,2> bounds{};
        std::copy_n(aabb.lowerBound.begin(), 3, bounds[0].begin());
        std::copy_n(aabb.upperBound.begin(), 3, bounds[1].begin());

        tmin = (bounds[mSign[0]][0] - mOrigin[0]) * mInvDirection[0];
        tmax = (bounds[1-mSign[0]][0] - mOrigin[0]) * mInvDirection[0];
        tymin = (bounds[mSign[1]][1] - mOrigin[1]) * mInvDirection[1];
        tymax = (bounds[1-mSign[1]][1] - mOrigin[1]) * mInvDirection[1];

        if ((tmin > tymax) || (tymin > tmax))
            return false;

        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;

        tzmin = (bounds[mSign[2]][2] -  mOrigin[2]) * mInvDirection[2];
        tzmax = (bounds[1-mSign[2]][2] -  mOrigin[2]) * mInvDirection[2];

        if ((tmin > tzmax) || (tzmin > tmax))
            return false;

        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;

        t = tmin;
        if (t < 0) {
            t = tmax;
            return false;
        }

        return true;
    }



inline
double deg2rad(const double &deg)
{ return deg * M_PI / 180; }

inline
double clamp(const double &lo, const double &hi, const double &v)
{ return std::max(lo, std::min(hi, v)); }

bool rayTriangleIntersect(
    const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
    double &t, double &u, double &v)
{
    //#ifdef MOLLER_TRUMBORE
    Vector3d v0v1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    Vector3d v0v2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

    // Cross product
    Vector3d pvec = {mDirection[1]*v0v2[2] - mDirection[2]*v0v2[1],
                     mDirection[2]*v0v2[0] - mDirection[0]*v0v2[2],
                     mDirection[0]*v0v2[1] - mDirection[1]*v0v2[0]};

    // Dot product
    double det = v0v1[0]*pvec[0] + v0v1[1]*pvec[1] + v0v1[2]*pvec[2];
    //#ifdef CULLING
        // if the determinant is negative the triangle is backfacing
        // if the determinant is close to 0, the ray misses the triangle
    if (fabs(det) < kEpsilon)
        return false;
    // #else
    //     // ray and triangle are parallel if det is close to 0
    //     if (fabs(det) < kEpsilon) return false;
    //#endif
    double invDet = 1 / det;

    // Subtraction
    Vector3d tvec = {mOrigin[0] - v0[0], mOrigin[1] - v0[1], mOrigin[2] - v0[2]};

    // Dot product x invDet
    u = (tvec[0]*pvec[0] + tvec[1]*pvec[1] + tvec[2]*pvec[2]) * invDet;
    //u = tvec.dotProduct(pvec) * invDet;

    if (u < 0 || u > 1)
        return false;

    // Cross product
    Vector3d qvec = {tvec[1]*v0v1[2] - tvec[2]*v0v1[1],
                     tvec[2]*v0v1[0] - tvec[0]*v0v1[2],
                     tvec[0]*v0v1[1] - tvec[1]*v0v1[0]};

    // Dot product x invDet
    v = (mDirection[0]*qvec[0] + mDirection[1]*qvec[1] + mDirection[2]*qvec[2]) * invDet;

    if (v < 0 || u + v > 1)
        return false;

    t = (v0v2[0]*qvec[0] + v0v2[1]*qvec[1] + v0v2[2]*qvec[2]) * invDet;

    if( t < 0 )
        return false;

    return true;
    // #else
    //     // compute plane's normal
    //     Vec3f v0v1 = v1 - v0;
    //     Vec3f v0v2 = v2 - v0;
    //     // no need to normalize
    //     Vec3f N = v0v1.crossProduct(v0v2);  //N
    //     double denom = N.dotProduct(N);

    //     // Step 1: finding P

    //     // check if ray and plane are parallel ?
    //     double NdotRayDirection = N.dotProduct(dir);

    //     if (fabs(NdotRayDirection) < kEpsilon)  //almost 0
    //         return false;  //they are parallel so they don't intersect !

    //     // compute d parameter using equation 2
    //     double d = -N.dotProduct(v0);

    //     // compute t (equation 3)
    //     t = -(N.dotProduct(orig) + d) / NdotRayDirection;

    //     // check if the triangle is in behind the ray
    //     if (t < 0) return false;  //the triangle is behind

    //     // compute the intersection point using equation 1
    //     Vec3f P = orig + t * dir;

    //     // Step 2: inside-outside test
    //     Vec3f C;  //vector perpendicular to triangle's plane

    //     // edge 0
    //     Vec3f edge0 = v1 - v0;
    //     Vec3f vp0 = P - v0;
    //     C = edge0.crossProduct(vp0);
    //     if (N.dotProduct(C) < 0) return false;  //P is on the right side

    //     // edge 1
    //     Vec3f edge1 = v2 - v1;
    //     Vec3f vp1 = P - v1;
    //     C = edge1.crossProduct(vp1);
    //     if ((u = N.dotProduct(C)) < 0)  return false;  //P is on the right side

    //     // edge 2
    //     Vec3f edge2 = v0 - v2;
    //     Vec3f vp2 = P - v2;
    //     C = edge2.crossProduct(vp2);
    //     if ((v = N.dotProduct(C)) < 0) return false;  //P is on the right side;

    //     u /= denom;
    //     v /= denom;

    //     return true;  //this ray hits the triangle
    // #endif
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
    AABB_tree(const TriangleMesh& TriangleMesh, double SkinThickness = 0.001, bool TouchIsOverlap=true) :
            aabb_base::Tree_base(3, SkinThickness, 16,  TouchIsOverlap)
    {

        // BaseTreeType::setBoxSize( {2, 2, 10} );
        // BaseTreeType::setPeriodicity( {false, false, false} );
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
        std::cout << "lower " << mLowerBound[0] << ", " << mLowerBound[1] << ", " << mLowerBound[2] << std::endl;
        std::cout << "upper " << mUpperBound[0] << ", " << mUpperBound[1] << ", " << mUpperBound[2] << std::endl;
        std::cout << "AABB_tree is build" << std::endl;
    }

    int count = 0;
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

    //std::vector<unsigned int> AABB_tree::query(unsigned int particle, const AABB_base& aabb_base)
    std::vector<unsigned int> query(Ray& rRay)
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
            count++;
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
            if (rRay.intersect(nodeAABB, t) )
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

    private:

    std::array<double,3> mLowerBound{};
    std::array<double,3> mUpperBound{};
    ///@}
};

///@}

#endif //AABB_tree_INCLUDE_H
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef AABB_PRIMITIVE_INCLUDE_H
#define AABB_PRIMITIVE_INCLUDE_H

/// External includes

/// External libraries
#include "aabb_tree/AABB_base.h"

/// Project includes

typedef aabb_base::AABB_base AABB_lohedges;

constexpr double kEpsilon = 1e-14;
constexpr double EpsilonBound = 1e-8;

///@name TIBRA Classes
///@{

/// Forward declaration
class AABB_primitive;

/**
 * @class  AABB_primitive
 * @author Manuel Messmer
 * @brief  Base class for aabb primitives.
*/
class AABB_primitive_base
{
public:
    ///@name Type Definitions
    ///@{
    typedef TriangleMesh::Vector3d Vector3d;

    ///@}
    ///@name Operations
    ///@{
    virtual bool intersect(const AABB_primitive &aabb) const = 0;

    virtual bool intersect(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                           double &t, double &u, double &v) const = 0;

    ///@}
};
///@} // end classes

/**
 * @class  AABB primitive
 * @author Manuel Messmer
 * @brief  AABB primitive to be used in AABB_tree.
*/
class AABB_primitive : public AABB_primitive_base, public AABB_lohedges
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef TriangleMesh::Vector3d Vector3d;
    typedef TriangleMesh::Vector3i Vector3i;

    ///@}
    ///@name Life cycle
    ///@{
    AABB_primitive(const Vector3d& rLowerBound, const Vector3d& rUpperBound) :
        AABB_lohedges(rLowerBound, rUpperBound)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if AABB intersects with aabb.
    ///@param aabb const AABB_primitive &aabb.
    ///@return bool.
    bool intersect(const AABB_primitive &aabb) const override {
        for (unsigned int i = 0; i < 3; ++i) {
            if (aabb.upperBound[i] <= lowerBound[i] || aabb.lowerBound[i] >= upperBound[i] ) {
                return false;
            }
        }
        return true;
    }

    ///@brief Returns true if AABB intersect with triangle.
    ///@brief This function uses the seperating axis theorem. https://dyn4j.org/2010/01/sat/
    ///@param v0 Vertex 1 of triangle.
    ///@param v1 Vertex 2 of triangle.
    ///@param v2 Vertex 3 of triangle.
    bool intersect(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                   double &t, double &u, double &v) const override {

        /// Get extent of aabb.
        Vector3d extent = { (upperBound[0] - lowerBound[0] )/2.0, (upperBound[1] - lowerBound[1] )/2.0, (upperBound[2] - lowerBound[2] )/2.0};


        /// Translate triangle to origin.
        Vector3d v0_orig = {v0[0] - centre[0], v0[1] - centre[1], v0[2] - centre[2]};
        Vector3d v1_orig = {v1[0] - centre[0], v1[1] - centre[1], v1[2] - centre[2]};
        Vector3d v2_orig = {v2[0] - centre[0], v2[1] - centre[1], v2[2] - centre[2]};

        // Compute the edge vectors of the triangle  (ABC). Line between vertices.
        Vector3d f0 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        Vector3d f1 = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
        Vector3d f2 = {v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};

        // Compute the face normals of the AABB, because the AABB
        // AABB is axis algined by definition.
        Vector3d u0 = {1.0, 0.0, 0.0};
        Vector3d u1 = {0.0, 1.0, 0.0};
        Vector3d u2 = {0.0, 0.0, 1.0};

        // There are a total of 13 axis to test!

        // First test (u0, u1, u2) vs. (f0, f1, f2). 9 tests in total.
        // u0 vs f0.
        Vector3d axis_u0_f0 = {u0[1]*f0[2] - u0[2]*f0[1], u0[2]*f0[0] - u0[0]*f0[2], u0[0]*f0[1] - u0[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f0) ){
            return false;
        }

        // u0 vs f1.
        Vector3d axis_u0_f1 = {u0[1]*f1[2] - u0[2]*f1[1], u0[2]*f1[0] - u0[0]*f1[2], u0[0]*f1[1] - u0[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f1) ){
            return false;
        }

        // u0 vs f2.
        Vector3d axis_u0_f2 = {u0[1]*f2[2] - u0[2]*f2[1], u0[2]*f2[0] - u0[0]*f2[2], u0[0]*f2[1] - u0[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u0_f2) ){
            return false;
        }

        // u1 vs f0.
        Vector3d axis_u1_f0 = {u1[1]*f0[2] - u1[2]*f0[1], u1[2]*f0[0] - u1[0]*f0[2], u1[0]*f0[1] - u1[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f0) ){
            return false;
        }

        // u1 vs f1.
        Vector3d axis_u1_f1 = {u1[1]*f1[2] - u1[2]*f1[1], u1[2]*f1[0] - u1[0]*f1[2], u1[0]*f1[1] - u1[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f1) ){
            return false;
        }

        // u1 vs f2.
        Vector3d axis_u1_f2 = {u1[1]*f2[2] - u1[2]*f2[1], u1[2]*f2[0] - u1[0]*f2[2], u1[0]*f2[1] - u1[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u1_f2) ){
            return false;
        }

        // u2 vs f0.
        Vector3d axis_u2_f0 = {u2[1]*f0[2] - u2[2]*f0[1], u2[2]*f0[0] - u2[0]*f0[2], u2[0]*f0[1] - u2[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f0) ){
            return false;
        }

        // u2 vs f1.
        Vector3d axis_u2_f1 = {u2[1]*f1[2] - u2[2]*f1[1], u2[2]*f1[0] - u2[0]*f1[2], u2[0]*f1[1] - u2[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f1) ){
            return false;
        }

        // u2 vs f2.
        Vector3d axis_u2_f2 = {u2[1]*f2[2] - u2[2]*f2[1], u2[2]*f2[0] - u2[0]*f2[2], u2[0]*f2[1] - u2[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, axis_u2_f2) ){
            return false;
        }

        // Test face normals of aabb. 3 Tests.
        // axis1: (1, 0, 0)
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u0) ){
            return false;
        }
        // axis2: (0, 1, 0)
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u1) ){
            return false;
        }

        // axis3 (0, 0, 1)
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, u2) ){
            return false;
        }

        // Test face normal of triangle
        Vector3d triangle_normal = {f0[1]*f1[2] - f0[2]*f1[1], f0[2]*f1[0] - f0[0]*f1[2], f0[0]*f1[1] - f0[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0_orig, v1_orig, v2_orig, extent, triangle_normal) ){
            return false;
        }

        // Return true of all checkes returned true.
        return true;
    }

private:

    ///@}
    ///@name Private Operations
    ///@{

    ///@brief Returns true if axis intersect.
    ///@details This function uses the seperating axis theorem.
    ///@param u0 Normal vector of aabb (1,0,0).
    ///@param u1 Normal vector of aabb (0,1,0).
    ///@param u2 Normal vector of aabb (0,0,1).
    ///@param v0 Vertex triangle 1. Triangle must be shifted to origin.
    ///@param v0 Vertex triangle 2. Triangle must be shifted to origin.
    ///@param v0 Vertex triangle 3. Triangle must be shifted to origin.
    ///@param extent Extent of aabb. Half of edge length in each direction.
    ///@param test_axis Axis to be tested.
    ///@return bool
    bool check_axis( const Vector3d &u0, const Vector3d &u1, const Vector3d &u2,
                     const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                     const Vector3d &extent, const Vector3d& test_axis ) const{

        // Project all 3 vertices of the triangle onto the test_axis.
        double pv0 = v0[0]*test_axis[0] + v0[1]*test_axis[1] + v0[2]*test_axis[2];
        double pv1 = v1[0]*test_axis[0] + v1[1]*test_axis[1] + v1[2]*test_axis[2];
        double pv2 = v2[0]*test_axis[0] + v2[1]*test_axis[1] + v2[2]*test_axis[2];

        // Project normals of aabb onto the test_axis.
        double pu0 = u0[0]*test_axis[0] + u0[1]*test_axis[1] + u0[2]*test_axis[2];
        double pu1 = u1[0]*test_axis[0] + u1[1]*test_axis[1] + u1[2]*test_axis[2];
        double pu2 = u2[0]*test_axis[0] + u2[1]*test_axis[1] + u2[2]*test_axis[2];

        // Project the AABB onto the seperating axis
        double r = extent[0] * std::abs(pu0) + extent[1] * std::abs(pu1) + extent[2] * std::abs(pu2);

        // Check if most extreme of the triangle points intersect r.
        if( std::max( -std::max({pv0, pv1, pv2}), std::min({pv0, pv1, pv2}) ) > r ){
            return false;
        }

        return true;
    }

    ///@}
};

/**
 * @class  Ray
 * @author Manuel Messmer
 * @brief  Ray to be used in AABB_tree. Direction of ray must be positive oriented.
 *         Check is omitted for better performance.
*/
class Ray_AABB_primitive : public AABB_primitive_base {

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
    Ray_AABB_primitive(const Vector3d& Origin, const Vector3d& Direction) :
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
    bool inside(const AABB_primitive &aabb) const
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
    ///@return bool
    bool intersect(const AABB_primitive &aabb) const override
    {
        // Return true if origin of Ray is inside the aabb
        if( inside(aabb) )
            return true;

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

        if (tmin < 0) {
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
        double &t, double &u, double &v) const override
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

        if (std::abs(det) < kEpsilon)
            return false;

        double invDet = 1 / det;

        // Substraction: mOrigin - v0
        Vector3d tvec = {mOrigin[0] - v0[0], mOrigin[1] - v0[1], mOrigin[2] - v0[2]};

        // Dot product x invDet: (tvec * pvec) * invDet
        u = (tvec[0]*pvec[0] + tvec[1]*pvec[1] + tvec[2]*pvec[2]) * invDet;

        if (u < -1e-10 || u > 1+1e-10)
            return false;

        // Cross product: tvec x v0v1
        Vector3d qvec = {tvec[1]*v0v1[2] - tvec[2]*v0v1[1],
                        tvec[2]*v0v1[0] - tvec[0]*v0v1[2],
                        tvec[0]*v0v1[1] - tvec[1]*v0v1[0]};

        // Dot product x invDet: (mDirection * qvec) * invDet
        v = (mDirection[0]*qvec[0] + mDirection[1]*qvec[1] + mDirection[2]*qvec[2]) * invDet;

        if (v < -1e-10 || u + v > 1+1e-10)
            return false;

        // Dot product: v0v2 * qvec
        t = (v0v2[0]*qvec[0] + v0v2[1]*qvec[1] + v0v2[2]*qvec[2]) * invDet;

        // Return false if ray intersects in negative direction.
        if( t < -1e-10 )
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


#endif // AABB_PRIMITIVE_INCLUDE_H
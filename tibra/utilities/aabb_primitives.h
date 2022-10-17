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


    bool intersect(const AABB_primitive &aabb) const override {
        for (unsigned int i = 0; i < 3; ++i) {
            if (aabb.upperBound[i] <= lowerBound[i] || aabb.lowerBound[i] >= upperBound[i]) {
                return false;
            }
        }

        return true;
    }

    bool check_axis( const Vector3d &u0, const Vector3d &u1, const Vector3d &u2,
                     const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                     const Vector3d &extent, const Vector3d& axis_u_f ) const{



        // // Testing axis: axis_u0_f0
        // // Project all 3 vertices of the triangle onto the Seperating axis
        // float p0 = Vector3.Dot(v0, axis_u0_f0);
        double pv0 = v0[0]*axis_u_f[0] + v0[1]*axis_u_f[1] + v0[2]*axis_u_f[2];
        // float p1 = Vector3.Dot(v1, axis_u0_f0);
        double pv1 = v1[0]*axis_u_f[0] + v1[1]*axis_u_f[1] + v1[2]*axis_u_f[2];
        // float p2 = Vector3.Dot(v2, axis_u0_f0);
        double pv2 = v2[0]*axis_u_f[0] + v2[1]*axis_u_f[1] + v2[2]*axis_u_f[2];

        double pu0 = u0[0]*axis_u_f[0] + u0[1]*axis_u_f[1] + u0[2]*axis_u_f[2];
        double pu1 = u1[0]*axis_u_f[0] + u1[1]*axis_u_f[1] + u1[2]*axis_u_f[2];
        double pu2 = u2[0]*axis_u_f[0] + u2[1]*axis_u_f[1] + u2[2]*axis_u_f[2];

        // Project the AABB onto the seperating axis
        // We don't care about the end points of the prjection
        // just the length of the half-size of the AABB
        // That is, we're only casting the extents onto the
        // seperating axis, not the AABB center. We don't
        // need to cast the center, because we know that the
        // aabb is at origin compared to the triangle!
        // float r = e.X * Math.Abs(Vector3.Dot(u0, axis_u0_f0)) +
        //             e.Y * Math.Abs(Vector3.Dot(u1, axis_u0_f0)) +
        //             e.Z * Math.Abs(Vector3.Dot(u2, axis_u0_f0));
        double r = extent[0] * std::abs(pu0) + extent[1] * std::abs(pu1) + extent[2] * std::abs(pu2);

        // // Now do the actual test, basically see if either of
        // // the most extreme of the triangle points intersects r
        // // You might need to write Min & Max functions that take 3 arguments
        if( std::max( -std::max({pv0, pv1, pv2}), std::min({pv0, pv1, pv2}) ) > r ){
            return false;
        }

        return true;
    }

    bool intersect(const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                   double &t, double &u, double &v) const override {


        // Vector3d e = aabb.Extents;
        Vector3d extent = { (upperBound[0] - lowerBound[0])/2.0, (upperBound[1] - lowerBound[1])/2.0, (upperBound[2] - lowerBound[2])/2.0};
        // // Translate the triangle as conceptually moving the AABB to origin
        // // This is the same as we did with the point in triangle test

        Vector3d v0c = {v0[0] - centre[0], v0[1] - centre[1], v0[2] - centre[2]};
        Vector3d v1c = {v1[0] - centre[0], v1[1] - centre[1], v1[2] - centre[2]};
        Vector3d v2c = {v2[0] - centre[0], v2[1] - centre[1], v2[2] - centre[2]};

        // // Compute the edge vectors of the triangle  (ABC)
        // // That is, get the lines between the points as vectors
        Vector3d f0 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
        Vector3d f1 = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
        Vector3d f2 = {v0[0] - v2[0], v0[1] - v2[1], v0[2] - v2[2]};

        // // Compute the face normals of the AABB, because the AABB
        // // is at center, and of course axis aligned, we know that
        // // it's normals are the X, Y and Z axis.

        Vector3d u0 = {1.0, 0.0, 0.0};
        Vector3d u1 = {0.0, 1.0, 0.0};
        Vector3d u2 = {0.0, 0.0, 1.0};

        // // There are a total of 13 axis to test!

        // // We first test against 9 axis, these axis are given by
        // // cross product combinations of the edges of the triangle
        // // and the edges of the AABB. You need to get an axis testing
        // // each of the 3 sides of the AABB against each of the 3 sides
        // // of the triangle. The result is 9 axis of seperation
        // // https://awwapp.com/b/umzoc8tiv/

        // // Compute the 9 axis
        // Vector3 axis_u0_f0 = Vector3.Cross(u0, f0);
        Vector3d axis_u0_f0 = {u0[1]*f0[2] - u0[2]*f0[1], u0[2]*f0[0] - u0[0]*f0[2], u0[0]*f0[1] - u0[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u0_f0) ){
            return false;
        }

        // Vector3 axis_u0_f1 = Vector3.Cross(u0, f1);
        Vector3d axis_u0_f1 = {u0[1]*f1[2] - u0[2]*f1[1], u0[2]*f1[0] - u0[0]*f1[2], u0[0]*f1[1] - u0[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u0_f1) ){
            return false;
        }

        // Vector3 axis_u0_f2 = Vector3.Cross(u0, f2);
        Vector3d axis_u0_f2 = {u0[1]*f2[2] - u0[2]*f2[1], u0[2]*f2[0] - u0[0]*f2[2], u0[0]*f2[1] - u0[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u0_f2) ){
            return false;
        }

        // Vector3 axis_u1_f0 = Vector3.Cross(u1, f0);
        Vector3d axis_u1_f0 = {u1[1]*f0[2] - u1[2]*f0[1], u1[2]*f0[0] - u1[0]*f0[2], u1[0]*f0[1] - u1[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u1_f0) ){
            return false;
        }

        // Vector3 axis_u1_f1 = Vector3.Cross(u1, f1);
        Vector3d axis_u1_f1 = {u1[1]*f1[2] - u1[2]*f1[1], u1[2]*f1[0] - u1[0]*f1[2], u1[0]*f1[1] - u1[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u1_f1) ){
            return false;
        }

        // Vector3 axis_u1_f2 = Vector3.Cross(u2, f2);
        Vector3d axis_u1_f2 = {u1[1]*f2[2] - u1[2]*f2[1], u1[2]*f2[0] - u1[0]*f2[2], u1[0]*f2[1] - u1[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u1_f2) ){
            return false;
        }

        // Vector3 axis_u2_f0 = Vector3.Cross(u2, f0);
        Vector3d axis_u2_f0 = {u2[1]*f0[2] - u2[2]*f0[1], u2[2]*f0[0] - u2[0]*f0[2], u2[0]*f0[1] - u2[1]*f0[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u2_f0) ){
            return false;
        }

        // Vector3 axis_u2_f1 = Vector3.Cross(u2, f1);
        Vector3d axis_u2_f1 = {u2[1]*f1[2] - u2[2]*f1[1], u2[2]*f1[0] - u2[0]*f1[2], u2[0]*f1[1] - u2[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u2_f1) ){
            return false;
        }

        // Vector3 axis_u2_f2 = Vector3.Cross(u2, f2);
        Vector3d axis_u2_f2 = {u2[1]*f2[2] - u2[2]*f2[1], u2[2]*f2[0] - u2[0]*f2[2], u2[0]*f2[1] - u2[1]*f2[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, axis_u2_f2) ){
            return false;
        }

        // // Next, we have 3 face normals from the AABB
        // // for these tests we are conceptually checking if the bounding box
        // // of the triangle intersects the bounding box of the AABB
        // // that is to say, the seperating axis for all tests are axis aligned:
        // // axis1: (1, 0, 0), axis2: (0, 1, 0), axis3 (0, 0, 1)
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, u0) ){
            return false;
        }

        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, u1) ){
            return false;
        }

        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, u2) ){
            return false;
        }
        // TODO: 3 SAT tests
        // // Do the SAT given the 3 primary axis of the AABB
        // // You already have vectors for this: u0, u1 & u2

        // // Finally, we have one last axis to test, the face normal of the triangle
        // // We can get the normal of the triangle by crossing the first two line segments
        Vector3d triangle_normal = {f0[1]*f1[2] - f0[2]*f1[1], f0[2]*f1[0] - f0[0]*f1[2], f0[0]*f1[1] - f0[1]*f1[0]};
        if( !check_axis(u0, u1, u2, v0c, v1c, v2c, extent, triangle_normal) ){
            return false;
        }
        // TODO: 1 SAT test

        // Passed testing for all 13 seperating axis that exist!
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


#endif // AABB_PRIMITIVE_INCLUDE_H
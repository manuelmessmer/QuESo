// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// Project includes
#include "embedding/ray_aabb_primitive.h"

// bool Ray_AABB_primitive::inside(const AABB_primitive &aabb) const {
//     if(    mOrigin[0] >= aabb.lowerBound[0]
//         || mOrigin[1] >= aabb.lowerBound[1]
//         || mOrigin[2] >= aabb.lowerBound[2]
//         || mOrigin[0] <= aabb.upperBound[0]
//         || mOrigin[1] <= aabb.upperBound[1]
//         || mOrigin[2] <= aabb.upperBound[2] ) {
//             return true;
//         }

//     return true;
// }

bool Ray_AABB_primitive::intersect(const AABB_primitive &aabb) const
{
    if( !mPositiveDir ){
        throw std::runtime_error( "Ray_AABB_primitive :: intersect :: Direction of Ray must be positive.");
    }

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    double lower_0 = aabb.lowerBound[0];
    double lower_1 = aabb.lowerBound[1];

    double upper_0 = aabb.upperBound[0];
    double upper_1 = aabb.upperBound[1];

    double origin_0 = mOrigin[0];
    double origin_1 = mOrigin[1];

    double inv_direction_0 = mInvDirection[0];
    double inv_direction_1 = mInvDirection[1];

    double lower_2 = aabb.lowerBound[2];
    double upper_2 = aabb.upperBound[2];
    double origin_2 = mOrigin[2];

    // Check if origin of lies inside aabb.
    if(    origin_0 >= lower_0
        && origin_1 >= lower_1
        && origin_2 >= lower_2
        && origin_0 <= upper_0
        && origin_1 <= upper_1
        && origin_2 <= upper_2 ) {
            return true;
        }

    tmin = (lower_0 - origin_0) * inv_direction_0;
    tymax = (upper_1 - origin_1) * inv_direction_1;
    if(tmin > tymax) {
        return false;
    }

    tmax = (upper_0 - origin_0) * inv_direction_0;
    tymin = (lower_1 - origin_1) * inv_direction_1;


    if( tymin > tmax )
        return false;

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

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

bool Ray_AABB_primitive::intersect( const Vector3d &v0, const Vector3d &v1, const Vector3d &v2,
                double &t, double &u, double &v, bool& BackFacing) const {

    // Substraction: v1-v0 and v2-v0
    Vector3d v0v1 = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]};
    Vector3d v0v2 = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]};

    // Cross product: mDirection x v0v2
    Vector3d pvec = {mDirection[1]*v0v2[2] - mDirection[2]*v0v2[1],
                        mDirection[2]*v0v2[0] - mDirection[0]*v0v2[2],
                        mDirection[0]*v0v2[1] - mDirection[1]*v0v2[0]};

    // Dot product: v0v1 * pvec
    double det = v0v1[0]*pvec[0] + v0v1[1]*pvec[1] + v0v1[2]*pvec[2];

    BackFacing = false;
    // If det is smaller than zero triangle is back facing.
    if( det < kEpsilon )
        BackFacing = true;

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
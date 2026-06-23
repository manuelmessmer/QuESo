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

#ifndef MATH_UTILITIES_HPP
#define MATH_UTILITIES_HPP

/// STL includes
#include <cmath>
//// Project includes
#include "queso/includes/define.hpp"

namespace queso {

namespace Math {

    ///@brief Simple Power function
    ///@param x Value
    ///@param p Order
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    inline double Pow( double x, std::size_t p){
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }

    ///@brief Dot product of two vectors.
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return double.
    inline double Dot(PointView rLhs, PointView rRhs) {
        return (rLhs[0]*rRhs[0] + rLhs[1]*rRhs[1] + rLhs[2]*rRhs[2]);
    }

    /// @brief Returns rLhs*rRhs (Elementwise)
    /// @param rLhs
    /// @param rRhs
    /// @return Vector3d
    inline Vector3d MultElementWise(PointView rLhs, PointView rRhs) {
        return {rLhs[0]*rRhs[0], rLhs[1]*rRhs[1], rLhs[2]*rRhs[2]};
    }

    ///@brief Cross product of two vectors.
    ///@param Vector3d rLHs
    ///@param Vector3d rRHs
    ///@return double.
    inline Vector3d Cross(PointView rLhs, PointView rRhs) {
        return { rLhs[1]*rRhs[2] - rLhs[2]*rRhs[1],
                 rLhs[2]*rRhs[0] - rLhs[0]*rRhs[2],
                 rLhs[0]*rRhs[1] - rLhs[1]*rRhs[0] };
    }

    ///@brief Norm of vector.
    ///@param Vector3d rLHs
    ///@return double.
    inline double Norm(PointView rLhs) {
        return std::sqrt( rLhs[0]*rLhs[0] + rLhs[1]*rLhs[1] + rLhs[2]*rLhs[2] );
    }

    /// @brief Returns max value of vector
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    inline T Max( const std::array<T, 3>& rVector ){
        return std::max<T>(std::max<T>(rVector[0], rVector[1]), rVector[2]);
    }

    /// @brief Returns min value of vector.
    /// @tparam T
    /// @param rVector
    /// @return T
    template<typename T>
    inline T Min( const std::array<T, 3>& rVector ){
        return std::min<T>(std::min<T>(rVector[0], rVector[1]), rVector[2]);
    }

} // End namespace math

inline Vector3d operator+(Vector3dView lhs, Vector3dView rhs)
{
    return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
}

inline Vector3d& operator+=(Vector3d& lhs, Vector3dView rhs)
{
    lhs[0] += rhs[0];
    lhs[1] += rhs[1];
    lhs[2] += rhs[2];
    return lhs;
}

inline Vector3i operator+(Vector3iView lhs, Vector3iView rhs)
{
    return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
}

inline Vector3d operator-(Vector3dView lhs, Vector3dView rhs)
{
    return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

inline Vector3d& operator-=(Vector3d& lhs, Vector3dView rhs)
{
    lhs[0] -= rhs[0];
    lhs[1] -= rhs[1];
    lhs[2] -= rhs[2];
    return lhs;
}

inline Vector3d operator*(double lhs, Vector3dView rhs)
{
    return {lhs * rhs[0], lhs * rhs[1], lhs * rhs[2]};
}

inline Vector3d operator*(Vector3dView lhs, double rhs)
{
    return {lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs};
}

inline Vector3d operator/(Vector3dView lhs, double rhs)
{
    return {lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs};
}

inline Vector3i operator-(Vector3iView lhs, Vector3iView rhs)
{
    return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

inline Vector3d& operator+=(Vector3d& lhs, double rhs)
{
    lhs[0] += rhs;
    lhs[1] += rhs;
    lhs[2] += rhs;
    return lhs;
}

inline Vector3d& operator*=(Vector3d& lhs, double rhs)
{
    lhs[0] *= rhs;
    lhs[1] *= rhs;
    lhs[2] *= rhs;
    return lhs;
}

inline Vector3d& operator/=(Vector3d& lhs, double rhs)
{
    lhs[0] /= rhs;
    lhs[1] /= rhs;
    lhs[2] /= rhs;
    return lhs;
}

inline Vector3i& operator*=(Vector3i& lhs, IndexType rhs)
{
    lhs[0] *= rhs;
    lhs[1] *= rhs;
    lhs[2] *= rhs;
    return lhs;
}

} // End namespace queso

#endif

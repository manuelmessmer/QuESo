// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef UTILITIES_H
#define UTILITIES_H

#include "containers/point_types.h"

namespace tibra {

namespace math {
    ///@brief Simple Power functions
    ///@param x Value
    ///@param p Order
    ///@details For gcc (without --ffast-math compiler flag) this is faster than std::pow().
    static double inline pow( double x, std::size_t p){
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }

    ///@brief Dot product of two vectors.
    ///@param PointType LHs
    ///@param PointType RHs
    ///@return double.
    static double inline dot( const PointType& rLhs, const PointType& rRhs) {
        return (rLhs[0]*rRhs[0] + rLhs[1]*rRhs[1] + rLhs[2]*rRhs[2]);
    }

} // End namespace math

} // End namespace tibra

#endif
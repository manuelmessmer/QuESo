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

#ifndef POLYNOMIAL_UTILITIES_INCLUDE_H
#define POLYNOMIAL_UTILITIES_INCLUDE_H

/// Project includes
#include "queso/includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/// Provides functions to evaluate Legendre Polynomials and their integrals.
namespace Polynomial {

    /// Simple power function
    inline double power( double x, std::size_t p) noexcept {
        double result = 1.0;
        while( p > 0UL ) {
            result = result * x;
            p -= 1;
        }
        return result;
    }

    /// Returns Legendre polynomial defined on (a,b)
    inline double f_x(double x, IndexType order, double a, double b) noexcept {
        switch(order) {
            case 0:
                return 1;
            case 1:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return tmp_x;
                }
            case 2:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/2.0*(3.0*power(tmp_x,2)-1.0);
                }
            case 3:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/2.0*(5.0*power(tmp_x,3) - 3.0*tmp_x);
                }
            case 4:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/8.0*(35.0*power(tmp_x,4)-30.0*power(tmp_x,2) +3.0);
                }
            case 5:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/8.0*(63.0*power(tmp_x,5)-70.0*power(tmp_x,3)+15.0*tmp_x);
                }
            case 6:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/16.0*(231.0*power(tmp_x,6)-315.0*power(tmp_x,4)+105.0*power(tmp_x,2)-5.0);
                }
            case 7:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/16.0*(429.0*power(tmp_x,7)-693.0*power(tmp_x,5)+315.0*power(tmp_x,3)-35.0*tmp_x);
                }
            case 8:
                {
                    double tmp_x = (2*x - a - b)/(b-a);
                    return 1.0/128.0*(6435.0*power(tmp_x,8) - 12012.0*power(tmp_x,6)+6930.0*power(tmp_x,4)-1260.0*power(tmp_x,2)+35.0);
                }
            default:
                return 0.0; // Do not throw. It will make the switch quite slow.
        };
    }


    /// Returns integral of Legendre polynomial defined on (a,b)
    inline double f_x_int(double x, IndexType order, double a, double b) noexcept {
        switch(order) {
            case 0:
                {
                    return x;
                }
            case 1:
                {
                    return -power((a + b - 2.0*x),2)/(4.0*(a - b));
                }
            case 2:
                {
                    return -x/2.0 - power((a + b - 2.0*x),3)/(4.0*power((a - b),2));
                }
            case 3:
                {
                    return (3.0*power( (a + b - 2.0*x),2) )/(8.0*(a - b)) - (5*power((a + b - 2.0*x),4))/(16*power((a - b),3));
                }
            case 4:
                {
                    return (3.0*x)/8.0 + (5.0*power((a + b - 2.0*x),3))/(8*power((a - b),2)) - (7.0*power((a + b - 2*x),5))/(16*power((a - b),4));
                }
            case 5:
                {
                    return (35*power((a + b - 2*x),4))/(32*power((a - b),3)) - (15*power((a + b - 2*x),2))/(32*(a - b)) - (21*power((a + b - 2*x),6))/(32*power((a - b),5));
                }
            case 6:
                {
                    return (63*power((a + b - 2*x),5))/(32*power((a - b),4)) - (35*power((a + b - 2*x),3))/(32*power((a - b),2)) - (5*x)/16 - (33*power((a + b - 2*x),7))/(32*power((a - b),6));
                }
            case 7:
                {
                    return (35.0*power( (a + b - 2*x),2))/(64*(a - b)) - (315*power((a + b - 2*x),4))/(128.0*power((a - b),3)) + (231.0*power((a + b - 2*x),6))/(64.0*power((a - b),5)) - (429.0*power((a + b - 2*x),8))/(256.0*power((a - b),7));
                }
            case 8:
                {
                    return (35.0*x)/128.0 + (105.0*power( (a + b - 2*x),3) )/(64.0*power( (a - b),2) ) - (693.0*power( (a + b - 2*x),5) )/(128.0*power((a - b),4) ) + (429.0*power((a + b - 2*x),7))/(64.0*power((a - b),6)) - (715.0*power( (a + b - 2*x),9) )/(256.0*power((a - b),8));
                }
            default:
                return 0.0; // Do not throw. It will make the switch quite slow.
        };
    }

} // End namespace Polynomial

} // End namespace queso

#endif // POLYNOMIAL_UTILITIES_INCLUDE_H

///@}


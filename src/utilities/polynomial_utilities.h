// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef POLYNOMIAL_UTILITIES_INCLUDE_H
#define POLYNOMIAL_UTILITIES_INCLUDE_H

// External includes
#include <variant>
#include <cmath>

namespace Polynomial {

    struct p0{};
    struct p1{};
    struct p2{};
    struct p3{};
    struct p4{};
    struct p5{};
    struct p6{};
    struct p7{};
    struct p8{};

    using Lp = std::variant<p0,p1,p2,p3,p4,p5,p6,p7,p8>;
    using Legendre = std::vector<Lp>;

    Legendre legendre{p0{}, p1{}, p2{}, p3{}, p4{}, p5{}, p6{}, p7{}, p8{}};

    struct F_x{
        // Shifted legendre polynomial (p=0)
        double inline operator()(const p0& p) const {
            return 1;
        }
        // Shifted legendre polynomial (p=1)
        double inline operator()(const p1& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return tmp_x;
        }
        // Shifted legendre polynomial (p=2)
        double inline operator()(const p2& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/2.0*(3.0*std::pow(tmp_x,2)-1.0);
        }
        // Shifted legendre polynomial (p=3)
        double inline operator()(const p3& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/2.0*(5.0*std::pow(tmp_x,3) - 3.0*tmp_x);
        }
        // Shifted legendre polynomial (p=4)
        double inline operator()(const p4& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/8.0*(35.0*std::pow(tmp_x,4)-30.0*std::pow(tmp_x,2) +3.0);
        }
        // Shifted legendre polynomial (p=5)
        double inline operator()(const p5& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/8.0*(63.0*std::pow(tmp_x,5)-70.0*std::pow(tmp_x,3)+15.0*tmp_x);
        }
        // Shifted legendre polynomial (p=6)
        double inline operator()(const p6& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/16.0*(231.0*std::pow(tmp_x,6)-315.0*std::pow(tmp_x,4)+105.0*std::pow(tmp_x,2)-5.0);
        }
        // Shifted legendre polynomial (p=7)
        double inline operator()(const p7& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/16.0*(429.0*std::pow(tmp_x,7)-693.0*std::pow(tmp_x,5)+315.0*std::pow(tmp_x,3)-35.0*tmp_x);
        }
        // Shifted legendre polynomial (p=8)
        double inline operator()(const p8& p) const {
            double tmp_x = (2*x - a - b)/(b-a);
            return 1.0/128.0*(6435.0*std::pow(tmp_x,8) - 12012.0*std::pow(tmp_x,6)+6930.0*std::pow(tmp_x,4)-1260.0*std::pow(tmp_x,2)+35.0);
        }
        double x{};
        double a{};
        double b{};
    };

double inline f_x_visit( const Lp& s, double x, double a, double b ) {
    return std::visit( F_x{x, a, b}, s );
}

// Caller for F_x
double inline f_x( double x, int order, double a, double b ){
    return f_x_visit( legendre[order], x, a, b);
}

struct F_x_int{
    // Integral of shifted legendre polynomial (p=0)
    double operator()(const p0& p) const {
        return x;
    }
    // Integral of shifted legendre polynomial (p=1)
    double operator()(const p1& p) const {
        return -std::pow((a + b - 2.0*x),2)/(4.0*(a - b));;
    }
    // Integral of shifted legendre polynomial (p=2)
    double operator()(const p2& p) const {
        return -x/2.0 - std::pow((a + b - 2.0*x),3)/(4.0*std::pow((a - b),2));
    }
    // Integral of shifted legendre polynomial (p=3)
    double operator()(const p3& p) const {
        return (3.0*std::pow( (a + b - 2.0*x),2) )/(8.0*(a - b)) - (5*std::pow((a + b - 2.0*x),4))/(16*std::pow((a - b),3));
    }
    // Integral of shifted legendre polynomial (p=4)
    double operator()(const p4& p) const {
        return (3.0*x)/8.0 + (5.0*std::pow((a + b - 2.0*x),3))/(8*std::pow((a - b),2)) - (7.0*std::pow((a + b - 2*x),5))/(16*std::pow((a - b),4));
    }
    // Integral of shifted legendre polynomial (p=5)
    double operator()(const p5& p) const {
        return (35*std::pow((a + b - 2*x),4))/(32*std::pow((a - b),3)) - (15*std::pow((a + b - 2*x),2))/(32*(a - b)) - (21*std::pow((a + b - 2*x),6))/(32*std::pow((a - b),5));
    }
    // Integral of shifted legendre polynomial (p=6)
    double operator()(const p6& p) const {
        return (63*std::pow((a + b - 2*x),5))/(32*std::pow((a - b),4)) - (35*std::pow((a + b - 2*x),3))/(32*std::pow((a - b),2)) - (5*x)/16 - (33*std::pow((a + b - 2*x),7))/(32*std::pow((a - b),6));
    }
    // Integral of shifted legendre polynomial (p=7)
    double operator()(const p7& p) const {
        return (35.0*std::pow( (a + b - 2*x),2))/(64*(a - b)) - (315*std::pow((a + b - 2*x),4))/(128.0*std::pow((a - b),3)) + (231.0*std::pow((a + b - 2*x),6))/(64.0*std::pow((a - b),5)) - (429.0*std::pow((a + b - 2*x),8))/(256.0*std::pow((a - b),7));
    }
    // Integral of shifted legendre polynomial (p=8)
    double operator()(const p8& p) const {
        return (35.0*x)/128.0 + (105.0*std::pow( (a + b - 2*x),3) )/(64.0*std::pow( (a - b),2) ) - (693.0*std::pow( (a + b - 2*x),5) )/(128.0*std::pow((a - b),4) ) + (429.0*std::pow((a + b - 2*x),7))/(64.0*std::pow((a - b),6)) - (715.0*std::pow( (a + b - 2*x),9) )/(256.0*std::pow((a - b),8));
    }
    double x{};
    double a{};
    double b{};
};

// Caller for F_x_int
double inline f_x_int_visit( const Lp& s, double x, double a, double b ){
    return std::visit( F_x_int{x, a, b}, s );
}

double inline f_x_int( double x, int order, double a, double b ){
    return f_x_int_visit( legendre[order], x, a, b);
}

} // End Namespace
#endif
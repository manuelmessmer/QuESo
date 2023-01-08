// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H
#define TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H

//// Projecy includes
#include "containers/integration_point.hpp"

namespace tibra {

class TriangleGaussLegendrePoints1
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static const PointArrayType& IntegrationPoints(){

        static const PointArrayType s_integration_points{{
            IntegrationPoint( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints2
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const PointArrayType& IntegrationPoints(){

        static const PointArrayType s_integration_points{{
            IntegrationPoint( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPoint( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPoint( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints3
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static const PointArrayType& IntegrationPoints(){

        static const PointArrayType s_integration_points{{
            IntegrationPoint( 1.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPoint( 3.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPoint( 1.00 / 5.00 , 3.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPoint( 1.00 / 3.00 , 1.00 / 3.00 , -27.00 / 96.00 )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints4
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const PointArrayType& IntegrationPoints(){
        const double wa = 0.054975871827661;
        const double wb = 0.1116907948390055;
        const double Na1 = 0.816847572980459;
        const double Nb1 = 0.108103018168070;
        const double Na2 = 0.091576213509771;
        const double Nb2 = 0.445948490915965;

        static const PointArrayType s_integration_points{{
            IntegrationPoint( Na2, Na2, wa ),
            IntegrationPoint( Na1, Na2, wa ),
            IntegrationPoint( Na2, Na1, wa ),
            IntegrationPoint( Nb2, Nb2, wb ),
            IntegrationPoint( Nb1, Nb2, wb ),
            IntegrationPoint( Nb2, Nb1, wb )
        }};

        return s_integration_points;
    }
};

} // End namespace tibra

#endif // TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H
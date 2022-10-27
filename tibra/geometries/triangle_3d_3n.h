// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_3D_3N
#define TRIANGLE_3D_3N

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <stdexcept>

// Project includes
#include "geometries/triangle_gauss_legendre_integration_points.h"

typedef std::array<std::array<double,2>,3> Matrix;
typedef std::vector<IntegrationPoint> IntegrationPointVectorType;
typedef std::unique_ptr<IntegrationPointVectorType> IntegrationPointVectorPtrType;
typedef std::array<double,3> PointType;
typedef std::size_t IndexType;
typedef std::size_t SizeType;

class Triangle3D3N
{
public:
    Triangle3D3N(PointType P1, PointType P2, PointType P3, PointType Normal) :
        mP1(P1), mP2(P2), mP3(P3), mNormalVector(Normal)
    {}

    Triangle3D3N(PointType P1, PointType P2, PointType P3) :
        mP1(P1), mP2(P2), mP3(P3)
    {}

    double Area()
    {
        const double a = std::sqrt( std::pow(mP1[0] - mP2[0], 2) + std::pow(mP1[1] - mP2[1], 2) + std::pow(mP1[2] - mP2[2], 2));
        const double b = std::sqrt( std::pow(mP2[0] - mP3[0], 2) + std::pow(mP2[1] - mP3[1], 2) + std::pow(mP2[2] - mP3[2], 2));
        const double c = std::sqrt( std::pow(mP3[0] - mP1[0], 2) + std::pow(mP3[1] - mP1[1], 2) + std::pow(mP3[2] - mP1[2], 2));

        const double s = (a+b+c) / 2.0;

        return std::sqrt(s*(s-a)*(s-b)*(s-c));
    }

    PointType& Normal(){
        // Check if normal is set..
        return mNormalVector;
    }

    PointType Center()
    {
        PointType tmp_point = {0.0, 0.0, 0.0};
        tmp_point[0] = 1.0/3.0 * (mP1[0] + mP2[0] + mP3[0]);
        tmp_point[1] = 1.0/3.0 * (mP1[1] + mP2[1] + mP3[1]);
        tmp_point[2] = 1.0/3.0 * (mP1[2] + mP2[2] + mP3[2]);

        return tmp_point;
    }

    PointType& P1(){
        return mP1;
    }

    PointType& P2(){
        return mP2;
    }

    PointType& P3(){
        return mP3;
    }

    Matrix& Jacobian(Matrix& rResults){
        rResults[0][0] = - mP1[0] + mP2[0];
        rResults[1][0] = - mP1[1] + mP2[1];
        rResults[2][0] = - mP1[2] + mP2[2];
        rResults[0][1] = - mP1[0] + mP3[0];
        rResults[1][1] = - mP1[1] + mP3[1];
        rResults[2][1] = - mP1[2] + mP3[2];

        return rResults;
    }

    static const std::vector<IntegrationPointVectorType>& AllIntegrationPoints()
    {
        static const std::vector<IntegrationPointVectorType> integration_points =
        {
            TriangleGaussLegendrePoints1::IntegrationPoints(),
            TriangleGaussLegendrePoints2::IntegrationPoints(),
            TriangleGaussLegendrePoints3::IntegrationPoints(),
            TriangleGaussLegendrePoints4::IntegrationPoints()
        };

        return integration_points;
    }

    // Todo: Add checl is method is larger than 3
    static const IntegrationPointVectorType& GetIntegrationPoints( IndexType Method ){

        return AllIntegrationPoints()[Method];
    }

    IntegrationPointVectorPtrType GetIntegrationPointsGlobal( IndexType Method ) {

        const auto& s_integration_points = GetIntegrationPoints(Method);
        const SizeType point_numbers = s_integration_points.size();

        IntegrationPointVectorPtrType p_global_integration_points
            = std::make_unique<IntegrationPointVectorType>(point_numbers);

        for( int i = 0; i < point_numbers; ++i){
            const double xx  = this->ShapeFunctionValue( 0, s_integration_points[i] ) * mP1[0] +
                               this->ShapeFunctionValue( 1, s_integration_points[i] ) * mP2[0] +
                               this->ShapeFunctionValue( 2, s_integration_points[i] ) * mP3[0] ;

            const double yy = this->ShapeFunctionValue( 0, s_integration_points[i] ) * mP1[1] +
                              this->ShapeFunctionValue( 1, s_integration_points[i] ) * mP2[1] +
                              this->ShapeFunctionValue( 2, s_integration_points[i] ) * mP3[1] ;

            const double zz = this->ShapeFunctionValue( 0, s_integration_points[i] ) * mP1[2] +
                              this->ShapeFunctionValue( 1, s_integration_points[i] ) * mP2[2] +
                              this->ShapeFunctionValue( 2, s_integration_points[i] ) * mP3[2] ;

            (*p_global_integration_points)[i] = IntegrationPoint(xx, yy, zz, s_integration_points[i].GetWeight() );
        }

        return std::move(p_global_integration_points);
    }

    double ShapeFunctionValue( IndexType ShapeFunctionIndex, const PointType& rPoint ) const {
        switch( ShapeFunctionIndex )
        {
        case 0:
            return( 1.0 -rPoint[0] - rPoint[1] );
        case 1:
            return( rPoint[0] );
        case 2:
            return( rPoint[1] );
        default:
            throw std::invalid_argument(" Triangle3D3N Wrong Index of Shape Function! ");
            break;
        }

        return 0;
    }



private:
    PointType mP1;
    PointType mP2;
    PointType mP3;

    PointType mNormalVector;
};

#endif
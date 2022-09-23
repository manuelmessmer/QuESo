// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// External includes
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <cmath>

// Project includes
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

typedef std::size_t SizeType;
typedef std::vector<std::array<double,2>> Ip1DVectorType;
typedef std::unique_ptr<Ip1DVectorType> Ip1DVectorPtrType;
typedef std::vector<std::vector<std::array<double, 2>>> Ip1DVectorVectorType;

// Public member functions
Ip1DVectorPtrType IntegrationPointFactory1D::GetGGQ( SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod Method ){
    if( Method == ReducedExact || Method == ReducedOrder1 || Method == ReducedOrder2){
        return GetGGQPoints(PolynomialDegree, NumberKnotSpans, Method);
    } else
    {
        throw std::invalid_argument("IntegrationPointFactory1D: Method not available");
    }
}

Ip1DVectorPtrType IntegrationPointFactory1D::GetGauss( SizeType PolynomialDegree, IntegrationMethod Method ){
    switch(Method)
    {
        case Gauss:
            return std::make_unique<Ip1DVectorType>(mGaussLegendrePoints[PolynomialDegree]);
        case ReducedGauss1:
            return std::make_unique<Ip1DVectorType>(mGaussLegendrePoints[PolynomialDegree-1]);
        case ReducedGauss2:
            return std::make_unique<Ip1DVectorType>(mGaussLegendrePoints[PolynomialDegree-2]);
        default:
            throw std::invalid_argument("IntegrationPointFactory1D: Method not available");
            break;
    }
}

const std::pair<SizeType, SizeType> IntegrationPointFactory1D::GetSpaceDimension(SizeType PolynomialDegre, IntegrationMethod Method ){
    switch(Method)
    {
        case ReducedExact:
            return {2*PolynomialDegre, PolynomialDegre-2};
        case ReducedOrder1:
            return {2*PolynomialDegre-1, PolynomialDegre-2};
        case ReducedOrder2:
            return {2*PolynomialDegre-2, PolynomialDegre-2};
        default:
            throw std::invalid_argument("IntegrationPointFactory1D: Method not available");
            break;
    }
}

Ip1DVectorPtrType IntegrationPointFactory1D::GetGGQPoints(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod Method){
    const double a = 0.0;
    const double b = 1.0;

    const auto dimension = GetSpaceDimension(PolynomialDegree, Method);

    const SizeType e = NumberKnotSpans;  // Number of elements
    const SizeType p = dimension.first;  // Degree
    const SizeType r = dimension.second; // Continuity

    const SizeType n = (p+1)*2 + (e-1)*(p-r) - p - 1; // Number of dofs
    const SizeType m = std::ceil(n/2.0);              // Number of quadrature points

    // Get correct base rule points
    Ip1DVectorVectorPtrType p_base_points{};
    if( p == 4 && r == 0 ){
        if( e % 2 == 0 ){
            p_base_points = S_4_0_base_even;
        }
        else {
            p_base_points = S_4_0_base_odd;
        }
    }
    else if(p == 6 && r == 2){
        if( e % 2 == 0 ){
            p_base_points = S_6_2_base_even;
        }
        else {
            p_base_points = S_6_2_base_odd;
        }
    }
    else if( p == 4 && r == 1 ){
        if( e % 2 == 0 ){
            int odd = m % 2;
            p_base_points = mBasePointsReduced2[PolynomialDegree-2][odd];
        }
        else {
            const SizeType odd = m % 2;
            if( odd ){
                p_base_points = S_4_1_base_odd;
            } else {
                // Requires sepcial rule!!
                p_base_points = S_4_1_base_even_2;
            }
        }
    }
    else {
        const SizeType odd = m % 2;
        switch(Method)
        {
            case ReducedExact:
                p_base_points = mBasePointsOptimal[PolynomialDegree-2][odd];
                break;
            case ReducedOrder1:
                p_base_points = mBasePointsReduced1[PolynomialDegree-2][odd];
                break;
            case ReducedOrder2:
                p_base_points = mBasePointsReduced2[PolynomialDegree-2][odd];
                break;
            default:
                throw std::invalid_argument("IntegrationPointFactory1D: Method not available1");
                break;
        }
    }

    SizeType m1 = (*p_base_points)[0].size(); // boundary nodes
    SizeType m2 = (*p_base_points)[1].size(); // internal nodes
    SizeType m3 = (*p_base_points)[2].size(); // center nodes

    Ip1DVectorType points(m);
    if( 2*m == n ){ // For odd number of nodes
        if( m > 2*m1 ){
            const Ip1DVectorVectorType& base_points =  *p_base_points;
            const SizeType z = std::ceil(0.5*m);
            std::copy_n(base_points[0].begin(), m1, points.begin());
            SizeType left = std::ceil(base_points[0][m1-1][0]);
            const SizeType right = std::ceil(base_points[1][m2-1][0]);
            SizeType ii = m1;
            while( ii < z ){
                std::copy_n(base_points[1].begin(), m2, points.begin()+ii);
                std::for_each(points.begin()+ii, points.begin()+ii+m2, [left](auto& rValue) { rValue[0] +=left;});
                left += right;
                ii += m2;
            }
            std::reverse_copy(points.begin(), points.begin()+z, points.end()- z );
            std::for_each(points.end()-z, points.end(), [e](auto& rValue) { rValue[0] = e - rValue[0];});
        }
        else {
            switch(Method)
            {
                case ReducedExact:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsOptimal[PolynomialDegree-2])[e-1]);
                case ReducedOrder1:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsReduced1[PolynomialDegree-2])[e-1]);
                case ReducedOrder2:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsReduced2[PolynomialDegree-2])[e-1]);
                default:
                    throw std::invalid_argument("IntegrationPointFactory1D: Method not available2");
                    break;
            }
        }
    }
    else { // For odd number of nodes
        if( m > 2*(m1+m3)-1 ){
            const Ip1DVectorVectorType& base_points =  *p_base_points;
            const SizeType z = std::ceil(0.5*m);
            std::copy_n(base_points[0].begin(), m1, points.begin());
            SizeType left = std::ceil(base_points[0].back()[0]);
            const SizeType right = std::ceil(base_points[1].back()[0]);
            SizeType ii = m1;
            while( ii < z ){
                std::copy_n(base_points[1].begin(), m2, points.begin()+ii);
                std::for_each(points.begin()+ii, points.begin()+ii+m2, [left](auto& rValue) { rValue[0] +=left;});
                left += right;
                ii += m2;
            }

            std::copy_n(base_points[2].begin(), m3, points.begin()+z-m3);
            std::for_each(points.begin()+z-m3, points.begin()+z, [e](auto& rValue) { rValue[0] +=0.5*e;});

            std::reverse_copy(points.begin(), points.begin()+z, points.end()- z );
            std::for_each(points.end()-z, points.end(), [e](auto& rValue) { rValue[0] = e - rValue[0];});
        }
        else {
            switch(Method)
            {
                case ReducedExact:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsOptimal[PolynomialDegree-2])[e-1]);
                case ReducedOrder1:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsReduced1[PolynomialDegree-2])[e-1]);
                case ReducedOrder2:
                    return std::make_unique<Ip1DVectorType>( (*mPrecomputedPointsReduced2[PolynomialDegree-2])[e-1]);
                default:
                    throw std::invalid_argument("IntegrationPointFactory1D: Method not available3");
                    break;
            }
        }
    }

    // Scale points to desired interval (a,b)
    const double h = (b-a) / e;
    std::for_each(points.begin(), points.end(), [a, h](auto& rValue) { rValue[0] = a + h*rValue[0];
                                                                        rValue[1] *= h; });

    return std::make_unique<Ip1DVectorType>(points);
}


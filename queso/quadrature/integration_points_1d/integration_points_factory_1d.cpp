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

//// STL includes
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <cmath>
//// Project includes
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

typedef IntegrationPointFactory1D::Ip1DVectorType Ip1DVectorType;
typedef IntegrationPointFactory1D::Ip1DVectorVectorType Ip1DVectorVectorType;
typedef IntegrationPointFactory1D::Ip1DVectorPtrType Ip1DVectorPtrType;

// Public member functions
Ip1DVectorPtrType IntegrationPointFactory1D::GetGGQ(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethodType Method ){
    const double a = 0.0;
    const double b = 1.0;

    const auto dimension = GetSpaceDimension(PolynomialDegree, Method);

    const SizeType e = NumberKnotSpans;  // Number of elements
    const SizeType p = dimension.first;  // Degree
    const SizeType r = dimension.second; // Continuity

    const SizeType n = (p+1)*2 + (e-1)*(p-r) - p - 1; // Number of dofs
    const SizeType m = static_cast<SizeType>( std::ceil(n/2.0) );              // Number of quadrature points

    // Get correct base rule points
    const Ip1DVectorVectorType& r_base_points = GetGGQBasePoints(PolynomialDegree, NumberKnotSpans, Method);

    SizeType m1 = r_base_points[0].size(); // boundary nodes
    SizeType m2 = r_base_points[1].size(); // internal nodes
    SizeType m3 = r_base_points[2].size(); // center nodes

    Ip1DVectorPtrType p_ggq_points = MakeUnique<Ip1DVectorType>(m);
    Ip1DVectorType& r_ggq_points = *p_ggq_points;
    if( 2*m == n ){ // For odd number of nodes
        if( m > 2*m1 ){
            const SizeType z = static_cast<SizeType>( std::ceil(0.5*m) );
            std::copy_n(r_base_points[0].begin(), m1, r_ggq_points.begin());
            SizeType left = static_cast<SizeType>( std::ceil(r_base_points[0][m1-1][0]) );
            const SizeType right = static_cast<SizeType>( std::ceil(r_base_points[1][m2-1][0]) );
            SizeType ii = m1;
            while( ii < z ){
                std::copy_n(r_base_points[1].begin(), m2, r_ggq_points.begin()+ii);
                std::for_each(r_ggq_points.begin()+ii, r_ggq_points.begin()+ii+m2, [left](auto& rValue) { rValue[0] +=left;});
                left += right;
                ii += m2;
            }
            std::reverse_copy(r_ggq_points.begin(), r_ggq_points.begin()+z, r_ggq_points.end()- z );
            std::for_each(r_ggq_points.end()-z, r_ggq_points.end(), [e](auto& rValue) { rValue[0] = e - rValue[0];});
        }
        else {
            switch(Method)
            {
                case IntegrationMethod::ggq_optimal:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsOptimal[PolynomialDegree-2][e-1]);
                case IntegrationMethod::ggq_reduced_1:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsReduced1[PolynomialDegree-2][e-1]);
                case IntegrationMethod::ggq_reduced_2:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsReduced2[PolynomialDegree-2][e-1]);
                default:
                    assert(false);
            }
        }
    }
    else { // For odd number of nodes
        if( m > 2*(m1+m3)-1 ){
            const SizeType z = static_cast<SizeType>( std::ceil(0.5*m) );
            std::copy_n(r_base_points[0].begin(), m1, r_ggq_points.begin());
            SizeType left = static_cast<SizeType>( std::ceil(r_base_points[0].back()[0]) );
            const SizeType right = static_cast<SizeType>( std::ceil(r_base_points[1].back()[0]) );
            SizeType ii = m1;
            while( ii < z ){
                std::copy_n(r_base_points[1].begin(), m2, r_ggq_points.begin()+ii);
                std::for_each(r_ggq_points.begin()+ii, r_ggq_points.begin()+ii+m2, [left](auto& rValue) { rValue[0] +=left;});
                left += right;
                ii += m2;
            }

            std::copy_n(r_base_points[2].begin(), m3, r_ggq_points.begin()+z-m3);
            std::for_each(r_ggq_points.begin()+z-m3, r_ggq_points.begin()+z, [e](auto& rValue) { rValue[0] +=0.5*e;});

            std::reverse_copy(r_ggq_points.begin(), r_ggq_points.begin()+z, r_ggq_points.end()- z );
            std::for_each(r_ggq_points.end()-z, r_ggq_points.end(), [e](auto& rValue) { rValue[0] = e - rValue[0];});
        }
        else {
            switch(Method)
            {
                case IntegrationMethod::ggq_optimal:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsOptimal[PolynomialDegree-2][e-1]);
                case IntegrationMethod::ggq_reduced_1:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsReduced1[PolynomialDegree-2][e-1]);
                case IntegrationMethod::ggq_reduced_2:
                    return MakeUnique<Ip1DVectorType>(mPrecomputedPointsReduced2[PolynomialDegree-2][e-1]);
                default:
                    assert(false);
            }
        }
    }

    // Scale points to desired interval (a,b)
    const double h = (b-a) / e;
    std::for_each(r_ggq_points.begin(), r_ggq_points.end(), [a, h](auto& rValue) { rValue[0] = a + h*rValue[0];
                                                                                   rValue[1] *= h; });
    return p_ggq_points;
}

const Ip1DVectorType& IntegrationPointFactory1D::GetGauss( SizeType PolynomialDegree, IntegrationMethodType Method ){
    switch(Method)
    {
        case IntegrationMethod::gauss:
            return mGaussLegendrePoints[PolynomialDegree];
        case IntegrationMethod::gauss_reduced_1:
            return mGaussLegendrePoints[PolynomialDegree-1];
        case IntegrationMethod::gauss_reduced_2:
            return mGaussLegendrePoints[PolynomialDegree-2];
        default:
            assert(false);
    }
}

const std::pair<SizeType, SizeType> IntegrationPointFactory1D::GetSpaceDimension(SizeType PolynomialDegre, IntegrationMethodType Method ){
    switch(Method)
    {
        case IntegrationMethod::ggq_optimal:
            return {2*PolynomialDegre, PolynomialDegre-2};
        case IntegrationMethod::ggq_reduced_1:
            return {2*PolynomialDegre-1, PolynomialDegre-2};
        case IntegrationMethod::ggq_reduced_2:
            return {2*PolynomialDegre-2, PolynomialDegre-2};
        default:
            assert(false);
    }
}

const Ip1DVectorVectorType& IntegrationPointFactory1D::GetGGQBasePoints(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethodType Method){

    const auto dimension = GetSpaceDimension(PolynomialDegree, Method);

    const SizeType e = NumberKnotSpans;  // Number of elements
    const SizeType p = dimension.first;  // Degree
    const SizeType r = dimension.second; // Continuity

    const SizeType n = (p+1)*2 + (e-1)*(p-r) - p - 1; // Number of dofs
    const SizeType m = static_cast<SizeType>( std::ceil(n/2.0) );              // Number of quadrature points

    // Get correct base rule points
    if( p == 4 && r == 0 ){
        if( e % 2 == 0 ){
            return S_4_0_base_even;
        }
        else {
            return S_4_0_base_odd;
        }
    }
    else if(p == 6 && r == 2){
        if( e % 2 == 0 ){
            return S_6_2_base_even;
        }
        else {
            return S_6_2_base_odd;
        }
    }
    else if( p == 4 && r == 1 ){
        if( e % 2 == 0 ){
            int odd = m % 2;
            return mBasePointsReduced2[PolynomialDegree-2][odd];
        }
        else {
            const SizeType odd = m % 2;
            if( odd ){
                return S_4_1_base_odd;
            } else {
                // Requires sepcial rule!!
                return S_4_1_base_even_2;
            }
        }
    }
    else {
        const SizeType odd = m % 2;
        switch(Method)
        {
            case IntegrationMethod::ggq_optimal:
                return mBasePointsOptimal[PolynomialDegree-2][odd];
            case IntegrationMethod::ggq_reduced_1:
                return mBasePointsReduced1[PolynomialDegree-2][odd];
            case IntegrationMethod::ggq_reduced_2:
                return mBasePointsReduced2[PolynomialDegree-2][odd];
            default:
                assert(false);
        }
    }
}
} // End namespace queso
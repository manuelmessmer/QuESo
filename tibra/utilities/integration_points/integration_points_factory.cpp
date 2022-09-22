// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#include "utilities/integration_points/integration_points_factory.h"

typedef std::size_t SizeType;
typedef std::vector<std::array<double,2>> Ip1DVectorType;
typedef std::unique_ptr<Ip1DVectorType> Ip1DVectorPtrType;
typedef std::vector<std::vector<std::array<double, 2>>> Ip1DVectorVectorType;

// Public member functions
Ip1DVectorPtrType IntegrationPointFactory::GetGGQ( int PolynomialDegree, int NumberKnotSpans, IntegrationMethod method ){
    if( method == ReducedExact || method == ReducedOrder1 || method == ReducedOrder2){
        return GetGGQPoints(PolynomialDegree, NumberKnotSpans, method);
    } else
    {
        throw std::invalid_argument("IntegrationPointFactory: Method not available");
    }

}

Ip1DVectorPtrType IntegrationPointFactory::GetGauss( int PolynomialDegree, IntegrationMethod method ){
    switch(method)
    {
        case Gauss:
            return std::make_unique<Ip1DVectorType>(IntegrationPoints::GaussLegendrePoints[PolynomialDegree-1]);
        case ReducedGauss1:
            return std::make_unique<Ip1DVectorType>(IntegrationPoints::GaussLegendrePoints[PolynomialDegree-2]);
        case ReducedGauss2:
            return std::make_unique<Ip1DVectorType>(IntegrationPoints::GaussLegendrePoints[PolynomialDegree-3]);
        default:
            throw std::invalid_argument("IntegrationPointFactory: Method not available");
            break;
    }
}

const std::pair<SizeType, SizeType> IntegrationPointFactory::GetSpaceDimension(SizeType PolynomialDegre, IntegrationMethod method ){

    switch(method)
    {
        case ReducedExact:
            return {2*PolynomialDegre, PolynomialDegre-2};
        case ReducedOrder1:
            return {2*PolynomialDegre-1, PolynomialDegre-2};
        case ReducedOrder2:
            return {2*PolynomialDegre-2, PolynomialDegre-2};
        default:
            throw std::invalid_argument("IntegrationPointFactory: Method not available");
            break;
    }
}

Ip1DVectorPtrType IntegrationPointFactory::GetGGQPoints(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod method){
    const double a = 0.0;
    const double b = 1.0;

    const auto dimension = GetSpaceDimension(PolynomialDegree, method);

    const SizeType e = NumberKnotSpans;
    const SizeType p = dimension.first;
    const SizeType r = dimension.second;

    const SizeType n = (p+1)*2 + (e-1)*(p-r) - p - 1;
    const SizeType m = std::ceil(n/2.0);


    Ip1DVectorType points(m);


    // Get correct base rule points
    std::vector<std::vector<std::array<double, 2>>> base_points{};
    if( p == 4 && r == 0 ){
        if( e % 2 == 0 ){
            base_points = *S_4_0_base_even;
        }
        else {
            base_points = *S_4_0_base_odd;
        }
    }
    else if(p == 6 && r == 2){
        if( e % 2 == 0 ){
            base_points = *S_6_2_base_even;
        }
        else {
            base_points = *S_6_2_base_odd;
        }
    }
    else if( p == 4 && r == 1 ){
        if( e % 2 == 0 ){
            int odd = m % 2;
            base_points = *base_points_reduced2[PolynomialDegree-2][odd];
        }
        else {
            int odd = m % 2;
            if( odd ){
                base_points = *S_4_1_base_odd;
            } else {
                // Requires sepcial rule!!
                base_points = *S_4_1_base_even_2;
            }
        }
    }
    else {
        int odd = m % 2;
        switch(method)
        {
            case ReducedExact:
                base_points = *base_points_optimal[PolynomialDegree-2][odd];
                break;
            case ReducedOrder1:
                base_points = *base_points_reduced1[PolynomialDegree-2][odd];
                break;
            case ReducedOrder2:
                base_points = *base_points_reduced2[PolynomialDegree-2][odd];
                break;
            default:
                throw std::invalid_argument("IntegrationPointFactory: Method not available1");
                break;
        }
    }

    SizeType m1 = base_points[0].size(); // boundary nodes
    SizeType m2 = base_points[1].size(); // internal nodes
    SizeType m3 = base_points[2].size(); // center nodes

    if( 2*m == n ){
        // TODO..
        if( m > 2*m1 ){
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
            switch(method)
            {
                case ReducedExact:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_optimal[PolynomialDegree-2])[e-1]);
                case ReducedOrder1:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_reduced1[PolynomialDegree-2])[e-1]);
                case ReducedOrder2:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_reduced2[PolynomialDegree-2])[e-1]);
                default:
                    throw std::invalid_argument("IntegrationPointFactory: Method not available2");
                    break;
            }
        }
    }
    else { // For odd number of nodes
        if( m > 2*(m1+m3)-1 ){
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
            switch(method)
            {
                case ReducedExact:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_optimal[PolynomialDegree-2])[e-1]);
                case ReducedOrder1:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_reduced1[PolynomialDegree-2])[e-1]);
                case ReducedOrder2:
                    return std::make_unique<Ip1DVectorType>( (*precomputed_points_reduced2[PolynomialDegree-2])[e-1]);
                default:
                    throw std::invalid_argument("IntegrationPointFactory: Method not available3");
                    break;
            }

        }
    }

    const double h = (b-a) / e;
    std::for_each(points.begin(), points.end(), [a, h](auto& rValue) { rValue[0] = a + h*rValue[0];
                                                                        rValue[1] *= h; });
    int count = 1;


    return std::make_unique<Ip1DVectorType>(points);
}


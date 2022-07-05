#ifndef INTEGRATION_POINTS_FACTORY_H
#define INTEGRATION_POINTS_FACTORY_H

#include "gauss_legendre_integration_points.h"
#include "knot_span_integration_points_p_2.h"
#include "knot_span_integration_points_p_3.h"
#include "knot_span_integration_points_p_4.h"

class IntegrationPointFactory {
public:
    // TODO: Change to enum!!
    enum IntegrationMethod {Gauss, ReducedGauss1, ReducedGauss2, ReducedExact, ReducedOrder1, ReducedOrder2};

    static const std::vector<std::array<double, 2>>& GetIntegrationPoints( int PolynomialDegree, int NumberKnotSpans, IntegrationMethod method ){
        switch(method)
        {
            case Gauss:
                return GetGaussLegendrePoints(PolynomialDegree);
            case ReducedGauss1:
                return GetGaussLegendrePoints(PolynomialDegree-1);
            case ReducedGauss2:
                return GetGaussLegendrePoints(PolynomialDegree-2);
            case ReducedExact:
                return AllMultiKnotSpanIntegrationPointsExact()[PolynomialDegree-2][NumberKnotSpans-1];
            case ReducedOrder1:
                return AllMultiKnotSpanIntegrationPointsReducedOrder1()[PolynomialDegree-2][NumberKnotSpans-1];
            case ReducedOrder2:
                return AllMultiKnotSpanIntegrationPointsReducedOrder2()[PolynomialDegree-2][NumberKnotSpans-1];
            default:
                throw std::runtime_error("IntegrationPointFactory: Method not available");
        }
    }

    static const std::vector<std::array<double, 2>>& GetIntegrationPoints( int PolynomialDegree, IntegrationMethod method ){
        switch(method)
        {
            case Gauss:
                return GetGaussLegendrePoints(PolynomialDegree);
            case ReducedGauss1:
                return GetGaussLegendrePoints(PolynomialDegree-1);
            case ReducedGauss2:
                return GetGaussLegendrePoints(PolynomialDegree-2);
            default:
                throw std::runtime_error("IntegrationPointFactory: Method not available");
        }
    }

private:

    static const std::vector<std::vector<std::vector<std::array<double, 2>>>>& AllMultiKnotSpanIntegrationPointsExact()
    {
        static const std::vector<std::vector<std::vector<std::array<double, 2>>>> integration_points =
        {
            IntegrationPoints::Points_S_4_0,
            IntegrationPoints::Points_S_6_1,
            IntegrationPoints::Points_S_8_2
        };

        return integration_points;
    }

    static const std::vector<std::vector<std::vector<std::array<double, 2>>>>& AllMultiKnotSpanIntegrationPointsReducedOrder1()
    {
        static const std::vector<std::vector<std::vector<std::array<double, 2>>>> integration_points =
        {
            IntegrationPoints::Points_S_3_0,
            IntegrationPoints::Points_S_5_1,
            IntegrationPoints::Points_S_7_2
        };

        return integration_points;
    }

    static const std::vector<std::vector<std::vector<std::array<double, 2>>>>& AllMultiKnotSpanIntegrationPointsReducedOrder2()
    {
        static const std::vector<std::vector<std::vector<std::array<double, 2>>>> integration_points =
        {
            IntegrationPoints::Points_S_2_0,
            IntegrationPoints::Points_S_4_1,
            IntegrationPoints::Points_S_6_2
        };

        return integration_points;
    }

    static const std::vector<std::array<double, 2>>& GetGaussLegendrePoints( int PolynomialDegree ){

        return IntegrationPoints::GaussLegendrePoints[PolynomialDegree-1];
    }
};
#endif // INTEGRATION_POINTS_FACTORY_H
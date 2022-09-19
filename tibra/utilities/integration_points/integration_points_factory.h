// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINTS_FACTORY_H
#define INTEGRATION_POINTS_FACTORY_H

// External includes
#include <cstddef>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <memory>
#include <utility>

#include <stdexcept>

// Project includes
#include "gauss_legendre_integration_points.h"
#include "knot_span_integration_points_p_2.h"
#include "knot_span_integration_points_p_3.h"
#include "knot_span_integration_points_p_4.h"

class IntegrationPointFactory {
public:
    // Typefef
    typedef std::size_t SizeType;
    typedef std::vector<std::array<double,2>> Ip1DVectorType;
    typedef std::unique_ptr<Ip1DVectorType> Ip1DVectorPtrType;

    typedef std::vector<std::vector<std::array<double, 2>>> Ip1DVectorVectorType;
    typedef std::shared_ptr<Ip1DVectorVectorType> Ip1DVectorVectorPtrType;
    // Enum Definition
    enum IntegrationMethod {Gauss, ReducedGauss1, ReducedGauss2, ReducedExact, ReducedOrder1, ReducedOrder2};

    IntegrationPointFactory() = default;
    IntegrationPointFactory(const IntegrationPointFactory &m) = delete;
    IntegrationPointFactory & operator= (const IntegrationPointFactory &) = delete;

    static Ip1DVectorPtrType GetGGQ( int PolynomialDegree, int NumberKnotSpans, IntegrationMethod method );

    static Ip1DVectorPtrType GetGauss( int PolynomialDegree, IntegrationMethod method );

private:

    //static const std::vector<Ip1DVectorVectorType>& PrecomputedPointsP2();
    //static void GetBasePoints(int m, int p, IntegrationMethod method);

    // static const std::vector<Ip1DVectorVectorType>& AllMultiKnotSpanIntegrationPointsReducedOrder1();

    // static const std::vector<Ip1DVectorVectorType>& AllMultiKnotSpanIntegrationPointsReducedOrder2();

    static const std::pair<SizeType, SizeType> GetSpaceDimension(SizeType PolynomialDegre, IntegrationMethod method );

    static Ip1DVectorPtrType GetGGQPoints(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod method);

    // Private Member variables
    // Precomputed points optimal
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> base_points_optimal;

    static const Ip1DVectorVectorPtrType S_4_0_base_even;

    static const Ip1DVectorVectorPtrType S_4_0_base_odd;

    static const Ip1DVectorVectorPtrType S_6_1_base_even;

    static const Ip1DVectorVectorPtrType S_6_1_base_odd;

    static const std::vector<Ip1DVectorVectorPtrType> precomputed_points_optimal;

    static const Ip1DVectorVectorPtrType S_4_0_precomputed;

    static const Ip1DVectorVectorPtrType S_6_1_precomputed;
};
#endif // INTEGRATION_POINTS_FACTORY_H
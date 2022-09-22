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
    // TODO: Add Gauss Legende points here!!
    // Base points (Optimal)
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> base_points_optimal;
    // Base optimal contains:
    // { { S_4_0_base_even, S_4_0_base_odd},
    //   { S_6_1_base_even, S_6_1_base_odd},
    //   { S_8_2_base_even, S_8_2_base_odd} }
    static const Ip1DVectorVectorPtrType S_4_0_base_even;
    static const Ip1DVectorVectorPtrType S_4_0_base_odd;
    static const Ip1DVectorVectorPtrType S_6_1_base_even;
    static const Ip1DVectorVectorPtrType S_6_1_base_odd;
    static const Ip1DVectorVectorPtrType S_8_2_base_even;
    static const Ip1DVectorVectorPtrType S_8_2_base_odd;
    // Precomputed points (Optimal)
    static const std::vector<Ip1DVectorVectorPtrType> precomputed_points_optimal;
    // Precomputed optimal contains:
    // {S_4_0_precomputed, S_6_1_precomputed, S_8_2_precomputed}
    static const Ip1DVectorVectorPtrType S_4_0_precomputed;
    static const Ip1DVectorVectorPtrType S_6_1_precomputed;
    static const Ip1DVectorVectorPtrType S_8_2_precomputed;

    // Base points (Reduced1)
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> base_points_reduced1;
    // Base optimal reduced1:
    // { { S_3_0_base_even, S_3_0_base_odd},
    //   { S_5_1_base_even, S_5_1_base_odd},
    //   { S_7_2_base_even, S_7_2_base_odd} }
    static const Ip1DVectorVectorPtrType S_3_0_base_even;
    static const Ip1DVectorVectorPtrType S_3_0_base_odd;
    static const Ip1DVectorVectorPtrType S_5_1_base_even;
    static const Ip1DVectorVectorPtrType S_5_1_base_odd;
    static const Ip1DVectorVectorPtrType S_7_2_base_even;
    static const Ip1DVectorVectorPtrType S_7_2_base_odd;
    // Precomputed points (Reduced1)
    static const std::vector<Ip1DVectorVectorPtrType> precomputed_points_reduced1;
    // Precomputed reduced contains:
    // {S_3_0_precomputed, S_5_1_precomputed, S_7_2_precomputed}
    static const Ip1DVectorVectorPtrType S_3_0_precomputed;
    static const Ip1DVectorVectorPtrType S_5_1_precomputed;
    static const Ip1DVectorVectorPtrType S_7_2_precomputed;

    // Base points (Reduced2)
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> base_points_reduced2;
    // Base optimal reduced2:
    // { { S_2_0_base_even, S_2_0_base_odd},
    //   { S_4_1_base_even, S_4_1_base_odd},
    //   { S_6_2_base_even, S_6_2_base_odd} }
    static const Ip1DVectorVectorPtrType S_2_0_base_even;
    static const Ip1DVectorVectorPtrType S_2_0_base_odd;
    static const Ip1DVectorVectorPtrType S_4_1_base_even;
    static const Ip1DVectorVectorPtrType S_4_1_base_even_2; // Required for odd element numbers, but even quadrature rules.
    static const Ip1DVectorVectorPtrType S_4_1_base_odd;
    static const Ip1DVectorVectorPtrType S_6_2_base_even;
    static const Ip1DVectorVectorPtrType S_6_2_base_odd;

    static const std::vector<Ip1DVectorVectorPtrType> precomputed_points_reduced2;
    // Precomputed reduced contains:
    // {S_2_0_precomputed, S_4_1_precomputed, S_6_2_precomputed}
    static const Ip1DVectorVectorPtrType S_2_0_precomputed;
    static const Ip1DVectorVectorPtrType S_4_1_precomputed;
    static const Ip1DVectorVectorPtrType S_6_2_precomputed;

};
#endif // INTEGRATION_POINTS_FACTORY_H
// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINTS_FACTORY_1D_H
#define INTEGRATION_POINTS_FACTORY_1D_H

// STL includes
#include <vector>
#include <array>
#include <memory>

namespace tibra {

///@name TIBRA Classes
///@{

////
/**
 * @class  IntegrationPointFactory1D
 * @author Manuel Messmer
 * @brief  Factory for 1D Integration points for single and multiple knot spans.
 * @details Available Quadrature rules:
 *          {Gauss, ReducedGauss1, ReducedGauss2, ReducedExact, ReducedOrder1, ReducedOrder2}
 * @todo Wrap TIBRA in namespace and put enum outside the class.
*/
class IntegrationPointFactory1D {
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t SizeType;
    typedef std::vector<std::array<double,2>> Ip1DVectorType;
    typedef std::unique_ptr<Ip1DVectorType> Ip1DVectorPtrType;

    typedef std::vector<std::vector<std::array<double, 2>>> Ip1DVectorVectorType;
    typedef std::shared_ptr<Ip1DVectorVectorType> Ip1DVectorVectorPtrType;

    ///@}
    ///@name  Enum's
    ///@{
    enum IntegrationMethod {Gauss, ReducedGauss1, ReducedGauss2, ReducedExact, ReducedOrder1, ReducedOrder2};

    ///@}
    ///@name Operations
    ///@{

    /// @brief Get Generalized Gaussian Quadrature (GGQ) 1D rules.
    /// @param PolynomialDegree
    /// @param NumberKnotSpans
    /// @param Method options - {ReducedExact, ReducedOrder1, ReducedOrder2}
    /// @return IntegrationPoint1DVectorPtrTypes. Points are defined on the interval (0,1).
    static Ip1DVectorPtrType GetGGQ( SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod Method );

    /// @brief Get standard 1D Gauss-Legendre quadrature rules.
    /// @param PolynomialDegree
    /// @param Method options - {Gauss, ReducedGauss1, ReducedGauss2
    /// @return IntegrationPoint1DVectorPtrTypes. Points are defined on the interval (0,1).
    static Ip1DVectorPtrType GetGauss( SizeType PolynomialDegree, IntegrationMethod Method );

    ///@}
private:

    /// @name Private operations
    ///@{

    /// @brief Get integration space dimension for GGQ.
    /// @param PolynomialDegre
    /// @param Method options - {ReducedExact, ReducedOrder1, ReducedOrder2}
    /// @return std::pair<SizeType, SizeType> {Degree, Continuity}
    static const std::pair<SizeType, SizeType> GetSpaceDimension(SizeType PolynomialDegre, IntegrationMethod Method );

    /// @brief   Constructs GGQ rules from precomputed rules.
    /// @details Algorithm taken from: R. Hiemstra et al. Optimal and reduced quadrature rules for tensor product and hierarchically
    ///          refined splines in isogeometric analysis. Comput. Methods Appl. Mech. Engrg. 316 (2017) 966–1004,
    ///          http://dx.doi.org/10.1016/j.cma.2016.10.049
    /// @param PolynomialDegree
    /// @param NumberKnotSpans
    /// @param Method options - {ReducedExact, ReducedOrder1, ReducedOrder2}
    /// @return InetgrationPoint1DVectorPtrType
    static Ip1DVectorPtrType GetGGQPoints(SizeType PolynomialDegree, SizeType NumberKnotSpans, IntegrationMethod Method);
    ///@}

    /// @name Private operations
    ///@{

    /// Standard Gauss Legendre points, index=[p-1]
    static const Ip1DVectorVectorType mGaussLegendrePoints;

    /// Base points (GGQ Optimal)
    /// mBasePointsOptimal contains:
    /// { { S_4_0_base_even, S_4_0_base_odd},
    ///   { S_6_1_base_even, S_6_1_base_odd},
    ///   { S_8_2_base_even, S_8_2_base_odd} }
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> mBasePointsOptimal;
    static const Ip1DVectorVectorPtrType S_4_0_base_even;
    static const Ip1DVectorVectorPtrType S_4_0_base_odd;
    static const Ip1DVectorVectorPtrType S_6_1_base_even;
    static const Ip1DVectorVectorPtrType S_6_1_base_odd;
    static const Ip1DVectorVectorPtrType S_8_2_base_even;
    static const Ip1DVectorVectorPtrType S_8_2_base_odd;
    /// Precomputed points (GGQ Optimal)
    /// mPrecomputedPointsOptimal:
    /// {S_4_0_precomputed, S_6_1_precomputed, S_8_2_precomputed}
    static const std::vector<Ip1DVectorVectorPtrType> mPrecomputedPointsOptimal;
    static const Ip1DVectorVectorPtrType S_4_0_precomputed;
    static const Ip1DVectorVectorPtrType S_6_1_precomputed;
    static const Ip1DVectorVectorPtrType S_8_2_precomputed;

    /// Base points (GGQ Reduced1)
    /// mBasePointsReduced1 contains:
    /// { { S_3_0_base_even, S_3_0_base_odd},
    ///   { S_5_1_base_even, S_5_1_base_odd},
    ///   { S_7_2_base_even, S_7_2_base_odd} }
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> mBasePointsReduced1;
    static const Ip1DVectorVectorPtrType S_3_0_base_even;
    static const Ip1DVectorVectorPtrType S_3_0_base_odd;
    static const Ip1DVectorVectorPtrType S_5_1_base_even;
    static const Ip1DVectorVectorPtrType S_5_1_base_odd;
    static const Ip1DVectorVectorPtrType S_7_2_base_even;
    static const Ip1DVectorVectorPtrType S_7_2_base_odd;
    /// Precomputed points (GGQ Reduced1)
    /// mPrecomputedPointsReduced1:
    /// {S_3_0_precomputed, S_5_1_precomputed, S_7_2_precomputed}
    static const std::vector<Ip1DVectorVectorPtrType> mPrecomputedPointsReduced1;
    static const Ip1DVectorVectorPtrType S_3_0_precomputed;
    static const Ip1DVectorVectorPtrType S_5_1_precomputed;
    static const Ip1DVectorVectorPtrType S_7_2_precomputed;

    /// Base points (GGQ Reduced2)
    /// mBasePointsReduced2:
    /// { { S_2_0_base_even, S_2_0_base_odd},
    ///   { S_4_1_base_even, S_4_1_base_odd},
    ///   { S_6_2_base_even, S_6_2_base_odd} }
    static const std::vector<std::array<Ip1DVectorVectorPtrType,2>> mBasePointsReduced2;
    static const Ip1DVectorVectorPtrType S_2_0_base_even;
    static const Ip1DVectorVectorPtrType S_2_0_base_odd;
    static const Ip1DVectorVectorPtrType S_4_1_base_even;
    /// Required for odd element numbers, but even quadrature rules.
    static const Ip1DVectorVectorPtrType S_4_1_base_even_2;
    static const Ip1DVectorVectorPtrType S_4_1_base_odd;
    static const Ip1DVectorVectorPtrType S_6_2_base_even;
    static const Ip1DVectorVectorPtrType S_6_2_base_odd;

    /// Precomputed points (GGQ Reduced2)
    /// mPrecomputedPointsReduced2:
    /// {S_2_0_precomputed, S_4_1_precomputed, S_6_2_precomputed}
    static const std::vector<Ip1DVectorVectorPtrType> mPrecomputedPointsReduced2;
    static const Ip1DVectorVectorPtrType S_2_0_precomputed;
    static const Ip1DVectorVectorPtrType S_4_1_precomputed;
    static const Ip1DVectorVectorPtrType S_6_2_precomputed;

    ///@}

}; // End IntegrationPointFactory1D class
///@} // End classes
} // End namespace tibra

#endif // INTEGRATION_POINTS_FACTORY_1D_H

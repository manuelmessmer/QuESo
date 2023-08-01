// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MOMENT_FITTING_UTILITIES_INCLUDE_H
#define MOMENT_FITTING_UTILITIES_INCLUDE_H

//// STL includes
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <array>
#include <variant>
//// Project includes
#include "embedding/octree.h"
#include "containers/element.hpp"
#include "containers/boundary_integration_point.hpp"
#include "utilities/parameters.h"

namespace tibra {

///@name TIBRA Classes
///@{

/**
 * @class  QuadratureTrimmedElement.
 * @author Manuel Messmer
 * @brief  Provides functions to create integration rules for trimmed elements.
**/
class QuadratureTrimmedElement{
public:
    ///@name Type Definition
    ///@{
    typedef Element::IntegrationPointVectorType IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef std::vector<BoundaryIntegrationPoint> BoundaryIPsVectorType;
    typedef Unique<BoundaryIPsVectorType> BoundaryIPsVectorPtrType;
    typedef boost::numeric::ublas::vector<double> VectorType;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Creates integration points for trimmed domain.
    ///@details 1. Distributes initial integration points uniformly in trimmed domain.
    ///         2. Computes constant terms of moment fitting equation.
    ///         3. Solves moment fitting equation in iterative point elimination algorithm.
    /// See: M. Meßmer et. al: Efficient CAD-integrated isogeometric analysis of trimmed solids,
    ///      Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
    ///@param rElement
    ///@param rParam
    static double AssembleIPs(Element& rElement, const Parameters& rParam);

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /// @brief Distributes point within trimmed domain using an octree. In each leaf node, Gauss points according to rIntegrationOrder are generated.
    ///        Only points inside the trimmed domain are considered.
    ///        Every time this function is called the otree is refined and more points are distributed.
    /// @param[out] rIntegrationPoint
    /// @param rOctree
    /// @param MinNumPoints Minimum Number of Points
    /// @param rIntegrationOrder Order of Gauss quadrature.
    static void DistributeIntegrationPoints(IntegrationPointVectorType& rIntegrationPoint, Octree<TrimmedDomainBase>& rOctree, SizeType MinNumPoints, const Vector3i& rIntegrationOrder);

    /// @brief Computes constant terms of moment fitting equation via volume integration points.
    /// @param[out] rConstantTerms
    /// @param pIntegrationPoints (Unique<T>)
    /// @param rElement
    /// @param rParam
    static void ComputeConstantTerms(VectorType& rConstantTerms, const IntegrationPointVectorPtrType& pIntegrationPoints, const Element& rElement, const Parameters& rParam);

    /// @brief Computes constant terms of moment fitting equation via boundary integration points. This functions uses the divergence theorem
    //         to transform the respective volume integrals to countour/surface integrals.
    /// @param[out] rConstantTerms
    /// @param pBoundaryIPs (Unique<T>)
    /// @param rElement
    /// @param rParam
    static void ComputeConstantTerms(VectorType& rConstantTerms, const BoundaryIPsVectorPtrType& pBoundaryIPs, const Element& rElement, const Parameters& rParam);

    /// @brief Set-Up and solve moment fitting equation. Solve the moment fitting equation for the weights of the integration points.
    ///        Computed weights are directly assigned to rIntegrationPoint.
    /// @param rConstantTerms
    /// @param[out] rIntegrationPoint
    /// @param rElement
    /// @param rParam
    /// @return double Residual
    static double MomentFitting(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, const Element& rElement, const Parameters& rParam);

    /// @brief Start point elimination algorihtm. Final quadrature rule is stored in rElement.
    /// @param rConstantTerms
    /// @param rIntegrationPoint
    /// @param rElement
    /// @param rParam
    /// @return double Residual
    static double PointElimination(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, Element& rElement, const Parameters& rParam);
}; // End Class


} // End namespace tibra

#endif // MOMENT_FITTING_UTILITIES_INCLUDE_H
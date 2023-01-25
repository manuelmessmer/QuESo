// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef SINGLE_ELEMENT_H
#define SINGLE_ELEMENT_H

//// STL includes
#include <vector>
#include <array>
//// Project includes
#include "containers/element.hpp"
#include "containers/integration_point.hpp"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace tibra {

///@name TIBRA Classes
///@{

////
/**
 * @class  SingleElement. Provides assembly opeartions for tensor-product quadrature rules of single non-trimmed element.
 * @author Manuel Messmer
 * @brief  Provides assembly for 3D quadrature rules.
 * @details Available Quadrature rules:
 *          {Gauss, Gauss_Reduced1, Gauss_Reduced2}
*/
class SingleElement {

public:
        ///@name Type Definitions
        ///@{
        typedef std::vector<IntegrationPoint> IntegrationPointType;

        ///@}
        ///@name Operations
        ///@{

        /// @brief Assemble tensor product quadrature rules.
        /// @param rElement
        /// @param rParam
        static void AssembleIPs(Element& rElement, const Parameters& rParam);

        /// @brief Assemble tensor product quadrature rules.
        /// @note This functions clears rIntegrationPoints.
        /// @param[out] rIntegrationPoints
        /// @param rLowerBoundParam LowerBound of element in parametric space.
        /// @param rUpperBoundParam LowerBound of element in parametric space.
        /// @param rOrder Order of quadrature rule.
        /// @param Method Integration method: Default - Gauss.
        static void AssembleIPs(IntegrationPointType& rIntegrationPoints, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam,
                                const Vector3i& rOrder, IntegrationMethodType Method = IntegrationMethod::Gauss );
        ///@}

}; // End Class SingleElement

///@}

} // End namespace tibra

#endif // SINGLE_ELEMENT_H
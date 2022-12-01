// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef SINGLE_ELEMENT_H
#define SINGLE_ELEMENT_H

//// STL includes
#include <vector>
#include <array>
//// Project includes
#include "containers/element.h"
#include "containers/integration_point.h"
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
 *          {Gauss, ReducedGauss1, ReducedGauss2}
*/
class SingleElement {

public:
        ///@name Type Definitions
        ///@{
        typedef std::vector<IntegrationPoint> IntegrationPointType;

        ///@}
        ///@name Operations
        ///@{

        /// @brief Assemble tensor product quadrature rules
        /// @param rElement
        /// @param rParam
        static void AssembleIPs(Element& rElement, const Parameters& rParam);
        ///@}

}; // End Class SingleElement

///@}

} // End namespace tibra

#endif // SINGLE_ELEMENT_H
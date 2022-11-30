// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef SINGLE_ELEMENT_H
#define SINGLE_ELEMENT_H

#include <vector>
#include <array>

// Project includes
#include "geometries/element.h"
#include "geometries/integration_point.h"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

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
        typedef size_t SizeType;
        typedef std::array<double,3> PointType;
        typedef std::array<int,3> IntArrayType;
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

#endif // SINGLE_ELEMENT_H
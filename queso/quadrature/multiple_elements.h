// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MULTI_KNOTSPAN_BOXES_UTILITIES_H
#define MULTI_KNOTSPAN_BOXES_UTILITIES_H

//// Project includes
#include "containers/element_container.hpp"
#include "includes/parameters.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  QuadratureMultipleElements. Provides assembly opeartions for tensor-product quadrature rules of multiple non-trimmed elements.
 * @author Manuel Messmer
 * @brief  Provides assembly for 3D quadrature rules.
 * @details Available Quadrature rules:
 *          {GGQ_Optimal, GGQ_Reduced1, GGQ_Reduced2}
*/
class QuadratureMultipleElements {

public:
    ///@name Type Defintitions
    ///@{
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Operations
    ///@{

    /// @brief Assemble tensor product quadrature rules.
    /// @param rElements
    /// @param rNumberOfElements
    /// @param rIntegrationOrder
    /// @param Method Integration method
    static void AssembleIPs(ElementContainer& rElements, const Vector3i& rNumberOfElements, const Vector3i& rIntegrationOrder, IntegrationMethodType Method);

private:

    ///@todo Refactor this
    static void AssignNumberNeighbours(ElementContainer::ElementVectorPtrType& rElements, IndexType direction);

    static ElementContainer::ElementPtrType NextElement(ElementContainer& rElements, std::size_t id, bool& found, int direction );

    static void StoreIntegrationPoints(ElementContainer::ElementVectorPtrType& rElements, std::array<int,3>& rNumberKnotspans,
            const Vector3i& rIntegrationOrder, IntegrationMethodType Method);

    static bool AllElementsVisited(ElementContainer& rElements);
}; // End Class QuadratureMultipleElements

///@}

} // End namespace queso

#endif
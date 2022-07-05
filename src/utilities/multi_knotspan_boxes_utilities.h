// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MULTI_KNOTSPAN_BOXES_UTILITIES_H
#define MULTI_KNOTSPAN_BOXES_UTILITIES_H

#include "geometries/element_container.h"
#include "utilities/parameters.h"

typedef std::size_t IndexType;
typedef std::size_t SizeType;

namespace MultiKnotspanBoxesUtilities {

    void ComputeIntegrationPoints(ElementContainer& rElements, const Parameters& rParameter);

    void AssignNumberNeighbours(ElementContainer::ElementVectorPtrType& rElements, IndexType direction, const Parameters& rParameters);

    ElementContainer::ElementPtrType NextElement(ElementContainer& rElements, std::size_t id, bool& found, int direction );

    void StoreIntegrationPoints(ElementContainer::ElementVectorPtrType& rElements, std::array<int,3>& rNumberKnotspans, const Parameters& rParameters);

    bool AllElementsVisited(ElementContainer& rElements);
} // End Namespace MultiKnotspanBoxesUtilities

#endif
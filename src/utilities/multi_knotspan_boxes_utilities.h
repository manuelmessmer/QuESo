// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MULTI_KNOTSPAN_BOXES_UTILITIES_H
#define MULTI_KNOTSPAN_BOXES_UTILITIES_H

// Project includes
#include "geometries/element_container.h"
#include "utilities/parameters.h"

class MultiKnotspanBoxesUtilities {

public:
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    static void CreateIntegrationPointsNonTrimmed(ElementContainer& rElements, const Parameters& rParameter);

private:
    static void AssignNumberNeighbours(ElementContainer::ElementVectorPtrType& rElements, IndexType direction, const Parameters& rParameters);

    static ElementContainer::ElementPtrType NextElement(ElementContainer& rElements, std::size_t id, bool& found, int direction );

    static void StoreIntegrationPoints(ElementContainer::ElementVectorPtrType& rElements, std::array<int,3>& rNumberKnotspans, const Parameters& rParameters);

    static bool AllElementsVisited(ElementContainer& rElements);
}; // End Class MultiKnotspanBoxesUtilities

#endif
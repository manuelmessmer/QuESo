// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MULTI_KNOTSPAN_UTILITIES_H
#define MULTI_KNOTSPAN_UTILITIES_H

#include "geometries/element_container.h"
#include "utilities/parameters.h"

typedef std::size_t IndexType;
typedef std::size_t SizeType;

namespace MultiKnotspanUtilities {

    void ComputeIntegrationPoints(ElementContainer& rElements, const Parameters& rParameter);

    void Store1DIntegrationPoints(ElementContainer::ElementVectorPtrType& rElements, IndexType direction, const Parameters& rParameters);

    void Assemble1DIntegrationPoints(ElementContainer& rElements, const Parameters& rParameter);

} // End Namespace MultiKnotspanUtilities

#endif
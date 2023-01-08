// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

//// Project includes
#include "define.hpp"

namespace tibra {

namespace MappingUtilities {

PointType FromGlobalToLocalSpace( const PointType& rGlobalCoord,
                                  const PointType& rLowerPoint,
                                  const PointType& rUpperPoint);

PointType FromLocalToGlobalSpace( const PointType& rLocalCoord,
                                  const PointType& rLowerPoint,
                                  const PointType& rUpperPoint);

} // End namespace MappingUtilities

} // End namespace tibra

#endif
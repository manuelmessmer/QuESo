// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

//// Project includes
#include "define.hpp"

namespace tibra {

namespace Mapping {

PointType GlobalToParam( const PointType& rGlobalCoord,
                         const PointType& rLowerPoint,
                         const PointType& rUpperPoint);

PointType ParamToGlobal( const PointType& rLocalCoord,
                         const PointType& rLowerPoint,
                         const PointType& rUpperPoint);

} // End namespace Mapping

} // End namespace tibra

#endif